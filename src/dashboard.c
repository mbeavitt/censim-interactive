#include "dashboard.h"
#include "config.h"
#include "hist.h"
#include <raylib.h>
#include "raygui.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Internal histograms are binned FINE (DASH_INTERNAL_BINS) over their generous
// fixed ranges; the dashboard then re-aggregates the *populated* window down to a
// consistent number of on-screen bars at draw time (see draw_hist). Keeping the
// internal resolution high means a narrow data region still has plenty of fine
// bins to resolve, so the displayed detail no longer collapses to "a couple bars"
// when a metric only occupies a sliver of its range. DASH_MAXBINS bounds the
// snapshot/aggregation buffers and must be >= DASH_INTERNAL_BINS.
#define DASH_INTERNAL_BINS 480
#define DASH_MAXBINS       512

typedef struct {
    int   nbins, log_scale;
    float min, max;
    long  counts[DASH_MAXBINS];
    long  total, maxcount, underflow, overflow;
} HistSnap;

// Compact axis-number formatting into a small rotating buffer (so several calls
// can appear in one DrawText line).
static const char *fmt(float v) {
    static char bufs[6][24];
    static int k = 0;
    char *o = bufs[k = (k + 1) % 6];
    float a = v < 0 ? -v : v;
    if (v == 0)            snprintf(o, 24, "0");
    else if (a >= 1e5f)    snprintf(o, 24, "%.0e", v);
    else if (a >= 100.0f)  snprintf(o, 24, "%.0f", v);
    else if (a >= 10.0f)   snprintf(o, 24, "%.0f", v);
    else if (a >= 1.0f)    snprintf(o, 24, "%.1f", v);
    else                   snprintf(o, 24, "%.2f", v);
    return o;
}

static void snap(HistSnap *s, const Histogram *h) {
    s->nbins = h->nbins < DASH_MAXBINS ? h->nbins : DASH_MAXBINS;
    s->log_scale = h->log_scale;
    s->min = h->min; s->max = h->max;
    long mx = 0;
    for (int i = 0; i < s->nbins; i++) {
        s->counts[i] = h->counts[i];
        if (h->counts[i] > mx) mx = h->counts[i];
    }
    s->total = h->total;
    s->maxcount = mx;
    s->underflow = h->underflow;
    s->overflow = h->overflow;
}

// x-value at the median of the distribution, or NAN if the median falls outside
// the histogram's range. Accounts for under/overflow counts.
static float snap_median(const HistSnap *s) {
    if (!s || s->total == 0) return NAN;
    long target = s->total / 2;
    long cum = s->underflow;
    if (cum >= target) return NAN;  // median below the plotted range
    for (int i = 0; i < s->nbins; i++) {
        cum += s->counts[i];
        if (cum >= target) {
            float lo, hi;
            if (s->log_scale) {
                float L = log10f(s->min), H = log10f(s->max);
                lo = powf(10.0f, L + (H - L) * i / s->nbins);
                hi = powf(10.0f, L + (H - L) * (i + 1) / s->nbins);
                return sqrtf(lo * hi);
            }
            lo = s->min + (s->max - s->min) * i / s->nbins;
            hi = s->min + (s->max - s->min) * (i + 1) / s->nbins;
            return 0.5f * (lo + hi);
        }
    }
    return NAN;  // median above the plotted range
}

// x-value at the left edge of bin i (i in [0, nbins]), log-aware.
static float bin_edge(const HistSnap *s, int i) {
    if (s->log_scale) {
        float L = log10f(s->min), H = log10f(s->max);
        return powf(10.0f, L + (H - L) * i / s->nbins);
    }
    return s->min + (s->max - s->min) * i / s->nbins;
}

// Fraction [0,1] of value v across an explicit [lo,hi] window (log-aware).
static float val_frac(float v, float lo, float hi, int log_scale) {
    float x = v;
    if (log_scale) {
        if (v <= 0 || lo <= 0) return -1.0f;
        x = log10f(v); lo = log10f(lo); hi = log10f(hi);
    }
    if (hi <= lo) return -1.0f;
    return (x - lo) / (hi - lo);
}

// Choose the displayed bin window [b0,b1]. With autoscale, fit to the span of
// *visible* bins at the current Y scale: under log-Y even single-count tail bins
// render, so include every populated bin; under linear-Y tiny bars vanish, so
// include only bins above ~0.5% of the peak (otherwise the window stretches out to
// an invisible tail). Without autoscale, show the full fixed range.
static void window_bins(const HistSnap *s, int autoscale, int log_y, int *b0, int *b1) {
    *b0 = 0; *b1 = s->nbins - 1;
    if (!autoscale) return;
    long thresh = 1;
    if (!log_y) { thresh = s->maxcount / 200; if (thresh < 1) thresh = 1; }  // ~0.5% of peak
    int lo = -1, hi = -1;
    for (int i = 0; i < s->nbins; i++) if (s->counts[i] >= thresh) { if (lo < 0) lo = i; hi = i; }
    if (lo < 0) return;  // no in-range data; keep full range
    *b0 = lo; *b1 = hi;
}

static const Color BG    = (Color){15, 20, 15, 255};    // near-black, matches single view
static const Color GRID  = (Color){0, 90, 40, 255};     // neutral chrome (panel/status/legend)
static const Color BAR   = (Color){0, 200, 110, 255};   // default accent (panel chrome)
static const Color REAL  = (Color){255, 70, 70, 255};   // A. thaliana ghost
static const Color UNIF  = (Color){150, 150, 160, 255}; // uniform-model ghost (unused for now)
static const Color MED   = (Color){255, 255, 255, 255}; // live median of the plotted data
static const Color FIT   = (Color){255, 235, 90, 255};  // power-law-with-cutoff fit (block size)

// Solve the n x n system A x = rhs by Gaussian elimination with partial pivoting.
// Used for the weighted least-squares polynomial fits.
static int solve_sys(const double *A, const double *rhs, double *out, int n) {
    double M[10][10], b[10];
    if (n > 10) return 0;
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) M[i][j] = A[i*n + j];
        b[i] = rhs[i];
    }
    for(int i=0; i<n; i++) {
        int piv = i;
        for(int j=i+1; j<n; j++) {
            if (fabs(M[j][i]) > fabs(M[piv][i])) piv = j;
        }
        if (fabs(M[piv][i]) < 1e-12) return 0;
        if (piv != i) {
            for(int j=i; j<n; j++) { double t = M[i][j]; M[i][j] = M[piv][j]; M[piv][j] = t; }
            double t = b[i]; b[i] = b[piv]; b[piv] = t;
        }
        for(int j=i+1; j<n; j++) {
            double f = M[j][i] / M[i][i];
            for(int k=i; k<n; k++) M[j][k] -= f * M[i][k];
            b[j] -= f * b[i];
        }
    }
    for(int i=n-1; i>=0; i--) {
        out[i] = b[i];
        for(int j=i+1; j<n; j++) out[i] -= M[i][j] * out[j];
        out[i] /= M[i][i];
    }
    return 1;
}

// Per-plot accent palette: a cohesive "cool scientific" run from mint -> cyan ->
// sky -> aqua -> lime, anchored by a warm amber on the COMPOSITE summary metric.
// Each plot's chrome (title, border, axis text) and bars derive from its accent,
// echoing the glowing green/cyan/amber traces of the single-view mission control.
static const Color ACCENT[6] = {
    {0,   230, 140, 255},  // unique repeats / kb  - mint green
    {0,   210, 255, 255},  // HORs / kb            - cyan
    {70,  170, 255, 255},  // HOR block size       - sky blue
    {0,   220, 200, 255},  // HOR block gap        - aqua / teal
    {150, 230, 90,  255},  // HOR similarity       - lime
    {255, 190, 70,  255},  // composite HOR metric - amber (warm anchor)
};

// Scale a colour's brightness by f (alpha preserved), clamping each channel.
static Color shade(Color c, float f) {
    float r = c.r * f, g = c.g * f, b = c.b * f;
    if (r > 255) r = 255; if (g > 255) g = 255; if (b > 255) b = 255;
    if (r < 0)   r = 0;   if (g < 0)   g = 0;   if (b < 0)   b = 0;
    return (Color){(unsigned char)r, (unsigned char)g, (unsigned char)b, c.a};
}

// Left edge value of bin i for a raw Histogram (the reference uses these).
static float edge_at(const Histogram *h, int i) {
    if (h->log_scale) {
        float L = log10f(h->min), H = log10f(h->max);
        return powf(10.0f, L + (H - L) * i / h->nbins);
    }
    return h->min + (h->max - h->min) * i / h->nbins;
}

// Overlay the real-data distribution as a density curve (red line) over the plot,
// normalised to its own peak, using the same x-window and y scale as the sim bars.
static void draw_ref_curve(const Histogram *ref, float px, float py, float pw, float ph,
                           float lo_v, float hi_v, int log_x, bool log_y) {
    if (!ref || ref->total == 0) return;
    float dmax = 0.0f, dmin = 0.0f; int seen = 0;
    for (int i = 0; i < ref->nbins; i++) {
        if (ref->counts[i] == 0) continue;
        float w = log_x ? (edge_at(ref, i + 1) - edge_at(ref, i)) : 1.0f;
        float d = (float)ref->counts[i] / w;
        if (!seen || d > dmax) dmax = d;
        if (!seen || d < dmin) dmin = d;
        seen = 1;
    }
    if (dmax <= 0.0f) return;
    if (dmin <= 0.0f) dmin = dmax;
    float lmn = log10f(dmin), lmx = log10f(dmax);
    if (lmx - lmn < 1e-6f) lmn = lmx - 1.0f;

    int have = 0; float x0 = 0, y0 = 0;
    for (int i = 0; i < ref->nbins; i++) {
        float c = log_x ? sqrtf(edge_at(ref, i) * edge_at(ref, i + 1))
                        : 0.5f * (edge_at(ref, i) + edge_at(ref, i + 1));
        float f = val_frac(c, lo_v, hi_v, log_x);
        if (f < 0.0f || f > 1.0f) { have = 0; continue; }  // outside window: break line
        float w = log_x ? (edge_at(ref, i + 1) - edge_at(ref, i)) : 1.0f;
        float d = (float)ref->counts[i] / w;
        float h01 = (d <= 0.0f) ? 0.0f
                  : (log_y ? (log10f(d) - lmn) / (lmx - lmn) : d / dmax);
        if (h01 < 0.0f) h01 = 0.0f; else if (h01 > 1.0f) h01 = 1.0f;
        float X = px + f * pw, Y = py + ph * (1.0f - h01);
        if (have) DrawLineEx((Vector2){x0, y0}, (Vector2){X, Y}, 2.0f, REAL);
        x0 = X; y0 = Y; have = 1;
    }
}

static void draw_hist(Rectangle b, const char *title, const HistSnap *s,
                      bool log_y, bool autoscale, const Histogram *ref,
                      int target_bars, Color accent, int fit_type, int poly_order,
                      char *out_fit_text) {
    if (out_fit_text) out_fit_text[0] = '\0';
    Color axc = shade(accent, 0.55f);  // dimmed accent for border + axis text
    DrawRectangleRec(b, BG);
    DrawRectangleLinesEx(b, 1, axc);
    DrawText(title, (int)b.x + 7, (int)b.y + 5, 11, shade(accent, 0.35f));  // glow
    DrawText(title, (int)b.x + 6, (int)b.y + 4, 11, accent);
    char nbuf[40];
    snprintf(nbuf, sizeof(nbuf), "n=%ld", s ? s->total : 0);
    DrawText(nbuf, (int)b.x + 6, (int)b.y + 17, 9, axc);

    // Warn if any data fell outside the binned range (clipped to under/overflow).
    if (s && (s->underflow > 0 || s->overflow > 0)) {
        char ob[48];
        snprintf(ob, sizeof(ob), "clipped <%ld >%ld", s->underflow, s->overflow);
        DrawText(ob, (int)(b.x + b.width) - MeasureText(ob, 9) - 6, (int)b.y + 17, 9, (Color){255,140,0,255});
    }

    // Inner plot area: leave a left margin for y numbers and a bottom strip for x.
    float px = b.x + 40, py = b.y + 30, pw = b.width - 48, ph = b.height - 46;

    if (!s || s->maxcount == 0) {
        DrawText("awaiting data", (int)(px + pw/2 - 40), (int)(py + ph/2), 10, axc);
        return;
    }

    // displayed x-window (autoscaled to visible bins at this Y scale) and value range
    int b0, b1; window_bins(s, autoscale, log_y, &b0, &b1);
    // Keep the real-data distribution in view too, so the sim-vs-real gap stays
    // visible: extend the window to cover the reference's populated value range.
    if (ref && ref->total > 0) {
        int rlo = -1, rhi = -1;
        for (int i = 0; i < ref->nbins; i++)
            if (ref->counts[i] > 0) { if (rlo < 0) rlo = i; rhi = i; }
        if (rlo >= 0) {
            float flo = val_frac(edge_at(ref, rlo),     s->min, s->max, s->log_scale);
            float fhi = val_frac(edge_at(ref, rhi + 1), s->min, s->max, s->log_scale);
            if (flo >= 0.0f) { int bb = (int)(flo * s->nbins);
                if (bb < 0) bb = 0; else if (bb >= s->nbins) bb = s->nbins - 1; if (bb < b0) b0 = bb; }
            if (fhi >= 0.0f) { int bb = (int)(fhi * s->nbins);
                if (bb < 0) bb = 0; else if (bb >= s->nbins) bb = s->nbins - 1; if (bb > b1) b1 = bb; }
        }
    }
    // Breathing room: pad the window by ~10% of its span (min 1 bin) each side.
    if (autoscale) {
        int pad = (b1 - b0 + 1) / 10;
        if (pad < 1) pad = 1;
        b0 -= pad; if (b0 < 0) b0 = 0;
        b1 += pad; if (b1 >= s->nbins) b1 = s->nbins - 1;
    }
    int avail = b1 - b0 + 1;          // fine bins inside the displayed window
    float lo_v = bin_edge(s, b0), hi_v = bin_edge(s, b1 + 1);

    // Adaptive bin sizer: re-aggregate the window's fine bins into a consistent
    // number of on-screen bars so every plot shows comparable resolution, no matter
    // how wide its populated region is. We can only merge fine bins, never split
    // them, so cap the bar count at what's available (a genuinely sparse metric
    // simply shows fewer bars rather than a stretched handful).
    int bars = target_bars; if (bars < 1) bars = 1; if (bars > avail) bars = avail;

    // Per-bar plotted value is a DENSITY (count / linear x-width). On a LOG x-axis
    // log-spaced bins widen toward the right, so raw counts pile into a misleading
    // hump; density recovers the true shape (e.g. block gap: high at small gaps,
    // decaying to rare large ones, as in the paper). Linear x has equal widths, so
    // density == count up to scale. Computing it per *display* bar keeps the y-axis
    // labels consistent with the bars actually drawn.
    static float dens[DASH_MAXBINS];  // single-threaded UI: reused each frame
    float dmax_v = 0.0f, dmin_v = 0.0f;
    int seen = 0;
    // Polynomial fit accumulators (weighted least squares of ln density)
    int poly_order_actual = (fit_type == 1) ? poly_order : (fit_type == 3 ? 1 : 0);
    int n_coeffs = poly_order_actual + 1;
    if (n_coeffs > 10) n_coeffs = 10;
    double M_cheb[100] = {0}, rhs_cheb[10] = {0}; int nfit=0;
    double u_min = (lo_v > 0.0f) ? log(lo_v) : -10.0;
    double u_max = (hi_v > 0.0f) ? log(hi_v) : 10.0;

    // Normal / Gamma fit accumulators
    double mu = 0, stddev = 0;
    long total_for_normal = 0;
    if ((fit_type == 2 || fit_type == 4) && s) {
        double sx = 0, sx2 = 0;
        for (int i = 0; i < s->nbins; i++) {
            if (s->counts[i] > 0) {
                double x = s->log_scale ? sqrt((double)bin_edge(s, i) * bin_edge(s, i+1)) 
                                        : 0.5 * (bin_edge(s, i) + bin_edge(s, i+1));
                sx += x * s->counts[i];
                sx2 += x * x * s->counts[i];
                total_for_normal += s->counts[i];
            }
        }
        if (total_for_normal > 1) {
            mu = sx / total_for_normal;
            double var = (sx2 - sx * sx / total_for_normal) / (total_for_normal - 1);
            if (var > 0) stddev = sqrt(var);
        }
    }

    for (int j = 0; j < bars; j++) {
        int fs = b0 + (int)((long)j       * avail / bars);
        int fe = b0 + (int)((long)(j + 1) * avail / bars);
        if (fe <= fs) fe = fs + 1;
        long c = 0;
        for (int i = fs; i < fe && i <= b1; i++) c += s->counts[i];
        float w = s->log_scale ? (bin_edge(s, fe) - bin_edge(s, fs)) : (float)(fe - fs);
        float v = (w > 0.0f) ? (float)c / w : 0.0f;
        dens[j] = v;
        if (v > 0.0f) {
            if (!seen || v > dmax_v) dmax_v = v;
            if (!seen || v < dmin_v) dmin_v = v;
            seen = 1;
        }
        if ((fit_type == 1 || fit_type == 3) && v > 0.0f && c > 0) {
            // Bar centre (geometric for log x). Each populated bar gets EQUAL weight:
            double xc = s->log_scale ? sqrt((double)bin_edge(s, fs) * bin_edge(s, fe))
                                     : 0.5 * (bin_edge(s, fs) + bin_edge(s, fe));
            if (xc > 0.0) {
                double u = log(xc), y = log((double)v);
                double z = (u_max > u_min) ? 2.0 * (u - u_min) / (u_max - u_min) - 1.0 : 0.0;
                double T[10];
                T[0] = 1.0;
                T[1] = z;
                for (int k = 2; k < n_coeffs; k++) T[k] = 2.0 * z * T[k-1] - T[k-2];
                for (int row = 0; row < n_coeffs; row++) {
                    for (int col = 0; col < n_coeffs; col++) {
                        M_cheb[row * n_coeffs + col] += T[row] * T[col];
                    }
                    rhs_cheb[row] += T[row] * y;
                }
                nfit++;
            }
        }
    }
    if (dmax_v <= 0.0f) dmax_v = 1.0f;
    if (dmin_v <= 0.0f) dmin_v = dmax_v;
    float lminv = log10f(dmin_v), lmaxv = log10f(dmax_v);
    if (lmaxv - lminv < 1e-6f) lminv = lmaxv - 1.0f;  // single distinct value

    Color bar_top = shade(accent, 1.25f);  // brighter cap for a touch of dimension
    for (int j = 0; j < bars; j++) {
        if (dens[j] <= 0.0f) continue;
        float h01 = log_y ? (log10f(dens[j]) - lminv) / (lmaxv - lminv) : dens[j] / dmax_v;
        if (h01 < 0.0f) h01 = 0.0f; else if (h01 > 1.0f) h01 = 1.0f;
        float hh = ph * h01;
        if (log_y && hh < 1.0f) hh = 1.0f;  // keep the smallest log bar visible
        // Derive both x-edges from the same formula so adjacent bars share an exact
        // integer boundary -- this is what keeps the inter-bar gaps perfectly even
        // (independent int-casts of position and width used to jitter them).
        int xL = (int)(px + pw * (float)j       / bars);
        int xR = (int)(px + pw * (float)(j + 1) / bars);
        int w = xR - xL - 1; if (w < 1) w = 1;  // uniform 1px gap
        int yT = (int)(py + ph - hh);
        DrawRectangle(xL, yT, w, (int)hh, accent);
        DrawRectangle(xL, yT, w, 1, bar_top);
    }

    // y-axis: top = max plotted value, bottom = 0 (lin) or min (log); scale tag.
    // For log-x plots the value is a density (count per unit x), else a raw count.
    DrawText(fmt(dmax_v), (int)b.x + 3, (int)py - 4, 9, axc);
    DrawText(log_y ? fmt(dmin_v) : "0", (int)b.x + 3, (int)(py + ph - 8), 9, axc);
    DrawText(log_y ? "log" : "lin", (int)b.x + 3, (int)(py + ph/2), 9, axc);
    DrawText(s->log_scale ? "dens" : "cnt", (int)b.x + 3, (int)(py + ph/2 + 11), 8, shade(accent, 0.3f));

    // x-axis: 3 ticks across the displayed window (geometric mid for log)
    float xmid = s->log_scale ? sqrtf(lo_v * hi_v) : 0.5f * (lo_v + hi_v);
    int ytick = (int)(py + ph + 2);
    DrawText(fmt(lo_v), (int)px, ytick, 9, axc);
    const char *m = fmt(xmid); DrawText(m, (int)(px + pw/2 - MeasureText(m, 9)/2), ytick, 9, axc);
    const char *x = fmt(hi_v); DrawText(x, (int)(px + pw - MeasureText(x, 9)), ytick, 9, axc);
    const char *xs = s->log_scale ? "log x" : "lin x";
    DrawText(xs, (int)(px + pw - MeasureText(xs, 9)), (int)py - 2, 9, shade(accent, 0.3f));

    // real-data distribution overlay (red density curve)
    draw_ref_curve(ref, px, py, pw, ph, lo_v, hi_v, s->log_scale, log_y);

    // Polynomial overlays (Chebyshev or Power law)
    if ((fit_type == 1 || fit_type == 3) && nfit >= n_coeffs) {
        double cf[10];
        if (solve_sys(M_cheb, rhs_cheb, cf, n_coeffs)) {
            double sse = 0.0;
            for (int j = 0; j < bars; j++) {
                if (dens[j] <= 0.0f) continue;
                int fs = b0 + (int)((long)j       * avail / bars);
                int fe = b0 + (int)((long)(j + 1) * avail / bars);
                if (fe <= fs) fe = fs + 1;
                double xc = s->log_scale ? sqrt((double)bin_edge(s, fs) * bin_edge(s, fe))
                                         : 0.5 * (bin_edge(s, fs) + bin_edge(s, fe));
                if (xc > 0.0) {
                    double u = log(xc);
                    double z = (u_max > u_min) ? 2.0 * (u - u_min) / (u_max - u_min) - 1.0 : 0.0;
                    double T[10]; T[0] = 1.0; T[1] = z;
                    for (int k = 2; k < n_coeffs; k++) T[k] = 2.0 * z * T[k-1] - T[k-2];
                    double log_d = 0;
                    for (int i = 0; i < n_coeffs; i++) log_d += cf[i] * T[i];
                    double d_pred = exp(log_d);
                    double diff = dens[j] - d_pred;
                    sse += diff * diff;
                }
            }
            if (out_fit_text) {
                int pos = 0;
                for (int i = 0; i < n_coeffs; i++) {
                    int n = snprintf(out_fit_text + pos, 128 - pos, "c%d=%.2g ", i, cf[i]);
                    if (n > 0) pos += n;
                    if (pos >= 127) break;
                }
                snprintf(out_fit_text + pos, 128 - pos, "sse=%.2g", sse);
            }
            int hv = 0; float xprev = 0, yprev = 0;
            for (int k = 0; k <= 64; k++) {
                float f = (float)k / 64.0f;
                float xv = s->log_scale
                         ? powf(10.0f, log10f(lo_v) + (log10f(hi_v) - log10f(lo_v)) * f)
                         : lo_v + (hi_v - lo_v) * f;
                if (xv <= 0.0) continue;
                double u = log(xv);
                double z = (u_max > u_min) ? 2.0 * (u - u_min) / (u_max - u_min) - 1.0 : 0.0;
                double T[10];
                T[0] = 1.0;
                T[1] = z;
                for (int i = 2; i < n_coeffs; i++) T[i] = 2.0 * z * T[i-1] - T[i-2];
                double log_d = 0;
                for (int i = 0; i < n_coeffs; i++) log_d += cf[i] * T[i];
                double d = exp(log_d);
                if (!(d > 0.0)) { hv = 0; continue; }
                float h01 = log_y ? (log10f((float)d) - lminv) / (lmaxv - lminv)
                                  : (float)d / dmax_v;
                if (h01 < 0.0f) h01 = 0.0f; else if (h01 > 1.0f) h01 = 1.0f;
                float X = px + f * pw, Y = py + ph * (1.0f - h01);
                if (hv) DrawLineEx((Vector2){xprev, yprev}, (Vector2){X, Y}, 2.0f, FIT);
                xprev = X; yprev = Y; hv = 1;
            }
            char fb[40]; snprintf(fb, sizeof(fb), (fit_type == 1) ? "Chebyshev %d order fit" : "Power law fit", poly_order_actual);
            DrawText(fb, (int)px + 3, (int)py + 4, 10, FIT);
        }
    } else if ((fit_type == 2 || fit_type == 4) && stddev > 0.0 && total_for_normal > 0 && mu > 0.0) {
        int hv = 0; float xprev = 0, yprev = 0;
        double fine_bin_w = (s->max - s->min) / s->nbins;
        double theta = (stddev * stddev) / mu;
        double k_shape = (mu * mu) / (stddev * stddev);
        for (int k = 0; k <= 64; k++) {
            float f = (float)k / 64.0f;
            float xv = s->log_scale
                     ? powf(10.0f, log10f(lo_v) + (log10f(hi_v) - log10f(lo_v)) * f)
                     : lo_v + (hi_v - lo_v) * f;
            double pdf_x = 0;
            if (fit_type == 2) {
                double z_norm = (xv - mu) / stddev;
                pdf_x = (1.0 / (stddev * sqrt(2.0 * 3.14159265358979323846))) * exp(-0.5 * z_norm * z_norm);
            } else if (fit_type == 4) {
                if (xv > 0.0) {
                    double log_pdf = (k_shape - 1.0) * log(xv) - (xv / theta) - k_shape * log(theta) - lgamma(k_shape);
                    pdf_x = exp(log_pdf);
                }
            }
            double d = total_for_normal * pdf_x * (s->log_scale ? 1.0 : fine_bin_w);
            if (!(d > 0.0)) { hv = 0; continue; }
            float h01 = log_y ? (log10f((float)d) - lminv) / (lmaxv - lminv)
                              : (float)d / dmax_v;
            if (h01 < 0.0f) h01 = 0.0f; else if (h01 > 1.0f) h01 = 1.0f;
            float X = px + f * pw, Y = py + ph * (1.0f - h01);
            if (hv) DrawLineEx((Vector2){xprev, yprev}, (Vector2){X, Y}, 2.0f, FIT);
            xprev = X; yprev = Y; hv = 1;
        }
        char fb[60];
        if (fit_type == 2) snprintf(fb, sizeof(fb), "Normal fit (mu=%.2g, std=%.2g)", mu, stddev);
        else snprintf(fb, sizeof(fb), "Gamma fit (k=%.2g, th=%.2g)", k_shape, theta);
        DrawText(fb, (int)px + 3, (int)py + 4, 10, FIT);
        
        double sse = 0.0;
        for (int j = 0; j < bars; j++) {
            if (dens[j] <= 0.0f) continue;
            int fs = b0 + (int)((long)j       * avail / bars);
            int fe = b0 + (int)((long)(j + 1) * avail / bars);
            if (fe <= fs) fe = fs + 1;
            double xc = s->log_scale ? sqrt((double)bin_edge(s, fs) * bin_edge(s, fe))
                                     : 0.5 * (bin_edge(s, fs) + bin_edge(s, fe));
            double pdf_x = 0;
            if (fit_type == 2) {
                double z_norm = (xc - mu) / stddev;
                pdf_x = (1.0 / (stddev * sqrt(2.0 * 3.14159265358979323846))) * exp(-0.5 * z_norm * z_norm);
            } else if (fit_type == 4) {
                if (xc > 0.0) {
                    double log_pdf = (k_shape - 1.0) * log(xc) - (xc / theta) - k_shape * log(theta) - lgamma(k_shape);
                    pdf_x = exp(log_pdf);
                }
            }
            double d_pred = total_for_normal * pdf_x * (s->log_scale ? 1.0 : fine_bin_w);
            double diff = dens[j] - d_pred;
            sse += diff * diff;
        }
        
        if (out_fit_text) {
            if (fit_type == 2) snprintf(out_fit_text, 128, "mu=%.2g, std=%.2g, sse=%.2g", mu, stddev, sse);
            else snprintf(out_fit_text, 128, "k=%.2g, theta=%.2g, sse=%.2g", k_shape, theta, sse);
        }
    }

    // live median of the plotted data (dashed bright marker)
    float med = snap_median(s);
    if (!isnan(med)) {
        float f = val_frac(med, lo_v, hi_v, s->log_scale);
        if (f >= 0 && f <= 1) {
            float xx = px + f * pw;
            for (int yy = (int)py; yy < (int)(py + ph); yy += 7)
                DrawLineEx((Vector2){xx, (float)yy}, (Vector2){xx, (float)yy + 4}, 2.5f, MED);
            char mb[32]; snprintf(mb, sizeof(mb), "med %s", fmt(med));
            DrawText(mb, (int)xx + 2, (int)(py + ph - 11), 9, MED);
        }
    }
}

void dashboard_init(Dashboard *d) {
    memset(d, 0, sizeof(*d));
    // Full paper baseline defaults.
    d->f_num_traj    = 300.0f;
    d->f_initial     = (float)DEFAULT_TRAJECTORY_INITIAL_SIZE; // 15000
    d->f_target_gens = 6000000.0f;
    d->f_collapse    = (float)DEFAULT_COLLAPSE_THRESHOLD;      // 300
    d->f_indel_rate  = DEFAULT_INDEL_RATE;
    d->f_snp_rate    = DEFAULT_SNP_RATE;
    d->f_indel_size  = DEFAULT_INDEL_SIZE_LAMBDA;
    d->f_size_ratio  = DEFAULT_DUP_DEL_SIZE_RATIO;
    d->elastic       = false;   // free drift by default
    d->f_elasticity  = 0.15f;
    d->f_nbins = 60.0f;
    d->f_cheb_order = 6.0f;
    d->f_sweep_steps = 10.0f;
    d->f_sweep_max_min = 5.0f;
    for (int i = 0; i < 6; i++) {
        d->plot_log_y[i] = false;
        d->fit_text[i][0] = '\0';
    }
    d->autoscale_x = false;
    // Per-plot log-Y defaults (left->right, top->bottom): lin lin log / log lin log
    d->plot_log_y[0] = false;  // unique/kb   lin
    d->plot_log_y[1] = false;  // HORs/kb     lin
    d->plot_log_y[2] = true;   // block size  log
    d->plot_log_y[3] = true;   // block gap   log
    d->plot_log_y[4] = false;  // similarity  lin
    d->plot_log_y[5] = true;   // composite   log

    // Real-data reference: try $CENSIM_REFERENCE, ./reference.dat, then next to the
    // executable. Absent is fine -- the overlay just won't draw.
    reference_init(&d->ref);
    d->show_ref = true;
    const char *env = getenv("CENSIM_REFERENCE");
    if (!(env && reference_load(&d->ref, env)) && !reference_load(&d->ref, "reference.dat")) {
        char path[1024];
        snprintf(path, sizeof(path), "%sreference.dat", GetApplicationDirectory());
        reference_load(&d->ref, path);
    }
}

void dashboard_free(Dashboard *d) {
    if (d->has_batch) { batch_free(&d->batch); d->has_batch = false; }
    reference_free(&d->ref);
}

static void launch(Dashboard *d) {
    if (d->has_batch) { batch_free(&d->batch); d->has_batch = false; }

    SimParams p;
    batch_default_params(&p, true);  // start from unbounded; add elastic pull below if enabled
    p.indel_rate        = d->f_indel_rate;
    p.snp_rate          = d->f_snp_rate;
    p.indel_size_lambda = d->f_indel_size;
    p.dup_del_size_ratio = d->f_size_ratio;
    p.dup_bias          = 1.0f / (1.0f + d->f_size_ratio);  // coupled freq: bigger dups => rarer (zero net drift)
    p.collapse_threshold = (int)d->f_collapse;
    p.target_size       = (int)d->f_initial;          // elastic pulls toward start size
    p.elasticity        = d->elastic ? d->f_elasticity : 0.0f;

    BatchConfig cfg;
    cfg.num_trajectories   = (int)d->f_num_traj;
    cfg.initial_size       = (int)d->f_initial;
    cfg.target_generations = (long)d->f_target_gens;
    cfg.seed_base          = (unsigned int)time(NULL);
    cfg.base_params        = p;
    cfg.nbins              = DASH_INTERNAL_BINS;  // fine internal bins; display re-aggregates

    batch_init(&d->batch, cfg, 0);
    batch_start(&d->batch);
    d->has_batch = true;
    d->started = true;
}

// Dump every batch histogram to a tidy CSV at full internal resolution. One row
// per (metric, bin); per-metric metadata (range, scale, totals, under/overflow)
// rides in '#'-prefixed comment lines so analysis tools (pandas comment='#') skip
// them. Edges are reported as the bin really was binned -- log-spaced for log-scale
// metrics -- so bin_lo/bin_center/bin_hi reconstruct the axis exactly. Returns the
// written path in `out` (caller-sized), or leaves out[0]=0 on failure.
static void export_histograms(Dashboard *d, const char *forced_path, char *out, size_t out_sz) {
    out[0] = 0;
    if (!d->has_batch) return;

    time_t now = time(NULL);
    struct tm tm; localtime_r(&now, &tm);
    char stamp[32];
    strftime(stamp, sizeof(stamp), "%Y%m%d_%H%M%S", &tm);

    char path[256];
    if (forced_path) {
        snprintf(path, sizeof(path), "%s", forced_path);
    } else {
        snprintf(path, sizeof(path), "censim_hist_%s.csv", stamp);
    }

    FILE *f = fopen(path, "w");
    if (!f) return;

    Batch *b = &d->batch;
    pthread_mutex_lock(&b->lock);

    struct { const char *name; const Histogram *h; } cols[] = {
        {"unique_per_kb", &b->h_unique_per_kb},
        {"hors_per_kb",   &b->h_hors_per_kb},
        {"block_size",    &b->h_block_size},
        {"block_gap",     &b->h_block_gap},
        {"similarity",    &b->h_similarity},
        {"diversity",     &b->h_diversity},
        {"composite",     &b->h_composite},
        {"collapse_gen",  &b->h_collapse_gen},
    };
    int ncols = (int)(sizeof(cols) / sizeof(cols[0]));

    fprintf(f, "# censim histogram export %s\n", stamp);
    fprintf(f, "# trajectories=%d survived=%d collapsed=%d\n",
            b->cfg.num_trajectories, b->survived, b->collapsed_count);
    for (int c = 0; c < ncols; c++) {
        const Histogram *h = cols[c].h;
        fprintf(f, "# %s: min=%g max=%g nbins=%d log_scale=%d "
                   "total=%ld underflow=%ld overflow=%ld\n",
                cols[c].name, h->min, h->max, h->nbins, h->log_scale,
                h->total, h->underflow, h->overflow);
    }

    fprintf(f, "metric,bin,bin_lo,bin_center,bin_hi,count\n");
    for (int c = 0; c < ncols; c++) {
        const Histogram *h = cols[c].h;
        for (int i = 0; i < h->nbins; i++)
            fprintf(f, "%s,%d,%g,%g,%g,%ld\n", cols[c].name, i,
                    hist_bin_lo(h, i), hist_bin_center(h, i),
                    hist_bin_lo(h, i + 1), h->counts[i]);
    }

    pthread_mutex_unlock(&b->lock);
    fclose(f);
    snprintf(out, out_sz, "%s", path);
}

// A labelled slider row; returns the (possibly updated) value via the pointer.
static void slider_row(float panel_x, float y, float w, const char *label,
                       const char *valtext, float *v, float lo, float hi) {
    DrawText(label, (int)panel_x + 12, (int)y + 2, 14, LIGHTGRAY);
    GuiSlider((Rectangle){panel_x + 130, y, w, 18}, NULL, valtext, v, lo, hi);
}

static void slider_sweep_row(float panel_x, float y, float w, const char *label,
                             const char *valtext, float *v, float lo, float hi, bool *sweep) {
    if (sweep) {
        GuiCheckBox((Rectangle){panel_x + 6, y + 2, 14, 14}, "", sweep);
        DrawText(label, (int)panel_x + 26, (int)y + 2, 14, LIGHTGRAY);
    } else {
        DrawText(label, (int)panel_x + 12, (int)y + 2, 14, LIGHTGRAY);
    }
    GuiSlider((Rectangle){panel_x + 130, y, w, 18}, NULL, valtext, v, lo, hi);
}

#include <sys/stat.h>

static void sweep_set_params(Dashboard *d) {
    int steps = (int)d->f_sweep_steps;
    if (steps < 1) steps = 1;
    if (d->sweep_enabled[0]) d->f_indel_rate = 0.0f + 3.0f * d->sweep_indices[0] / steps;
    if (d->sweep_enabled[1]) d->f_snp_rate = 0.0f + 1.0f * d->sweep_indices[1] / steps;
    if (d->sweep_enabled[2]) d->f_indel_size = 1.0f + 99.0f * d->sweep_indices[2] / steps;
    if (d->sweep_enabled[3]) {
        float r = -3.0f + 6.0f * d->sweep_indices[3] / steps;
        d->f_size_ratio = powf(10.0f, r);
    }
    if (d->sweep_enabled[4]) d->f_collapse = 0.0f + 5000.0f * d->sweep_indices[4] / steps;
    if (d->sweep_enabled[5]) d->f_elasticity = 0.0f + 1.0f * d->sweep_indices[5] / steps;
}

static void sweep_start(Dashboard *d) {
    int steps = (int)d->f_sweep_steps;
    if (steps < 1) steps = 1;
    d->sweep_total_runs = 1;
    for (int i = 0; i < 6; i++) {
        d->sweep_indices[i] = 0;
        d->sweep_counts[i] = d->sweep_enabled[i] ? (steps + 1) : 1;
        d->sweep_total_runs *= d->sweep_counts[i];
    }
    d->sweep_current_run = 1;
    d->sweep_running = true;
    
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    snprintf(d->sweep_dir, sizeof(d->sweep_dir), "sweep_%04d%02d%02d_%02d%02d%02d",
             tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min, tm->tm_sec);
    mkdir(d->sweep_dir, 0755);
    
    char path[512];
    snprintf(path, sizeof(path), "%s/histograms", d->sweep_dir);
    mkdir(path, 0755);
    
    snprintf(path, sizeof(path), "%s/summary.csv", d->sweep_dir);
    d->sweep_summary_file = fopen(path, "w");
    if (d->sweep_summary_file) {
        FILE *f = (FILE*)d->sweep_summary_file;
        fprintf(f, "run,indel_rate,snp_rate,indel_size,size_ratio,collapse,elasticity");
        fprintf(f, ",med_uniq,med_hors,med_size,med_gap,med_sim,med_comp");
        fprintf(f, ",fit_uniq,fit_hors,fit_size,fit_gap,fit_sim,fit_comp\n");
        fflush(f);
    }
    
    sweep_set_params(d);
    d->sweep_batch_start_time = GetTime();
    launch(d);
}

static void sweep_log_and_advance(Dashboard *d, const HistSnap *su, const HistSnap *sh, const HistSnap *bs, const HistSnap *bg, const HistSnap *si, const HistSnap *cp) {
    if (d->sweep_summary_file) {
        FILE *f = (FILE*)d->sweep_summary_file;
        fprintf(f, "%d,%.3f,%.3f,%.3f,%.3f,%.0f,%.3f",
                d->sweep_current_run, d->f_indel_rate, d->f_snp_rate, d->f_indel_size,
                d->f_size_ratio, d->f_collapse, d->f_elasticity);
        fprintf(f, ",%.2g,%.2g,%.2g,%.2g,%.2g,%.2g",
                snap_median(su), snap_median(sh), snap_median(bs), snap_median(bg), snap_median(si), snap_median(cp));
        for (int i=0; i<6; i++) {
            fprintf(f, ",\"%s\"", d->fit_text[i]);
        }
        fprintf(f, "\n");
        fflush(f);
    }
    
    char hpath[512], out[256];
    snprintf(hpath, sizeof(hpath), "%s/histograms/run_%04d.csv", d->sweep_dir, d->sweep_current_run);
    export_histograms(d, hpath, out, sizeof(out));

    d->sweep_current_run++;
    bool carry = true;
    for (int i = 0; i < 6; i++) {
        if (carry) {
            d->sweep_indices[i]++;
            if (d->sweep_indices[i] >= d->sweep_counts[i]) {
                d->sweep_indices[i] = 0;
                carry = true;
            } else {
                carry = false;
            }
        }
    }
    
    if (carry) {
        d->sweep_running = false;
        if (d->sweep_summary_file) {
            fclose((FILE*)d->sweep_summary_file);
            d->sweep_summary_file = NULL;
        }
    } else {
        sweep_set_params(d);
        d->sweep_batch_start_time = GetTime();
        launch(d);
    }
}

void dashboard_update_draw(Dashboard *d, int screen_w, int screen_h, int panel_w) {
    int panel_x = screen_w - panel_w;
    int plot_area_w = panel_x - 20;

    // ----- snapshot batch state under lock -----
    HistSnap su = {0}, sh = {0}, bs = {0}, bg = {0}, si = {0}, dv = {0}, cp = {0}, cg = {0};
    int completed = 0, survived = 0, collapsed = 0, total = 0, workers = 0, workers_live = 0;
    bool running = false, complete = false, stopping = false;
    if (d->has_batch) {
        Batch *b = &d->batch;
        total = b->cfg.num_trajectories;
        workers = b->num_workers;
        pthread_mutex_lock(&b->lock);
        completed = b->completed; survived = b->survived; collapsed = b->collapsed_count;
        workers_live = b->workers_running;
        snap(&su, &b->h_unique_per_kb); snap(&sh, &b->h_hors_per_kb);
        snap(&bs, &b->h_block_size);    snap(&bg, &b->h_block_gap);
        snap(&si, &b->h_similarity);    snap(&dv, &b->h_diversity);
        snap(&cp, &b->h_composite);     snap(&cg, &b->h_collapse_gen);
        pthread_mutex_unlock(&b->lock);
        // Reap once every worker has exited (natural finish OR after a stop request).
        // batch_join is instant here since the threads are already gone.
        if (b->running && workers_live == 0) batch_join(b);
        running   = b->running && !b->stop_requested;
        stopping  = b->stop_requested && workers_live > 0;
        complete  = !b->running && !b->stop_requested && completed >= total;
    }

    // ----- plot grid (left) -----
    int cols = 3, rows = 2;
    int margin = 8;
    int top = 64;  // leave room for status bar + mode toggle
    int gw = (plot_area_w - margin * (cols + 1)) / cols;
    int gh = (screen_h - top - margin * (rows + 1)) / rows;

    Rectangle cell[6];
    for (int i = 0; i < 6; i++) {
        int c = i % cols, r = i / cols;
        cell[i] = (Rectangle){ margin + c * (gw + margin),
                               top + margin + r * (gh + margin), gw, gh };
    }

    // With elastic bounds on, trajectories rarely collapse, so the "(survivors)"
    // qualifier on the per-array plots is just noise -- drop it in that mode.
    struct { const char *title; HistSnap *s; const char *metric; } plots[6] = {
        {d->elastic ? "UNIQUE REPEATS / kb" : "UNIQUE REPEATS / kb (survivors)", &su, "unique_per_kb"},
        {d->elastic ? "HORs / kb"           : "HORs / kb (survivors)",           &sh, "hors_per_kb"},
        {"HOR BLOCK SIZE",                  &bs, "block_size"},
        {"HOR BLOCK GAP",                   &bg, "block_gap"},
        {"HOR SIMILARITY",                  &si, "similarity"},
        {"COMPOSITE HOR METRIC",            &cp, "composite"},
    };
    Vector2 mouse = GetMousePosition();
    int target_bars = (int)d->f_nbins;
    for (int i = 0; i < 6; i++) {
        Color accent = ACCENT[i];
        const Histogram *ref = d->show_ref ? reference_get(&d->ref, plots[i].metric) : NULL;
        int fit_type = 0;
        if (i == 0) fit_type = 2; // Normal
        else if (i == 1) fit_type = 4; // Gamma
        else if (i == 2 || i == 3 || i == 5) fit_type = 1; // Chebyshev 6th
        else if (i == 4) fit_type = 3; // Power law
        draw_hist(cell[i], plots[i].title, plots[i].s, d->plot_log_y[i], d->autoscale_x,
                  ref, target_bars, accent, fit_type, (int)d->f_cheb_order,
                  d->fit_text[i]);
        // per-plot clickable Y-scale toggle (top-right of the cell)
        Rectangle tg = { cell[i].x + cell[i].width - 52, cell[i].y + 3, 46, 15 };
        bool hov = CheckCollisionPointRec(mouse, tg);
        DrawRectangleRec(tg, hov ? shade(accent, 0.3f) : (Color){25,35,25,255});
        DrawRectangleLinesEx(tg, 1, shade(accent, 0.55f));
        DrawText(d->plot_log_y[i] ? "Y:log" : "Y:lin", (int)tg.x + 5, (int)tg.y + 3, 9, accent);
        if (hov && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) d->plot_log_y[i] = !d->plot_log_y[i];
    }

    // ----- status bar -----
    DrawRectangle(0, top - 2, plot_area_w + 4, 2, GRID);
    
    if (d->sweep_running) {
        char params[256] = {0};
        int pos = snprintf(params, sizeof(params), "Testing:");
        if (d->sweep_enabled[0]) pos += snprintf(params+pos, sizeof(params)-pos, " INDEL rate=%.2f", d->f_indel_rate);
        if (d->sweep_enabled[1]) pos += snprintf(params+pos, sizeof(params)-pos, " SNP rate=%.2f", d->f_snp_rate);
        if (d->sweep_enabled[2]) pos += snprintf(params+pos, sizeof(params)-pos, " INDEL size=%.1f", d->f_indel_size);
        if (d->sweep_enabled[3]) pos += snprintf(params+pos, sizeof(params)-pos, " Dup/del ratio=%.2fx", d->f_size_ratio);
        if (d->sweep_enabled[4]) pos += snprintf(params+pos, sizeof(params)-pos, " Collapse=%.0f", d->f_collapse);
        if (d->sweep_enabled[5]) pos += snprintf(params+pos, sizeof(params)-pos, " Elasticity=%.2f", d->f_elasticity);
        DrawText(params, 8, 8, 14, (Color){150, 200, 255, 255});
    }

    char status[256];
    if (d->sweep_running) {
        float elapsed = GetTime() - d->sweep_batch_start_time;
        float rate = total ? 100.0f * collapsed / total : 0.0f;
        snprintf(status, sizeof(status),
                 "SWEEP [%d/%d]  |  %d/%d done   survived %d   collapsed %d (%.0f%%)   |  %d workers  |  run time: %.0fs",
                 d->sweep_current_run, d->sweep_total_runs, completed, total, survived, collapsed, rate, workers, elapsed);
    } else if (d->has_batch) {
        float rate = total ? 100.0f * collapsed / total : 0.0f;
        const char *state = stopping ? "STOPPING..." : running ? "RUNNING"
                          : complete ? "COMPLETE" : "STOPPED";
        snprintf(status, sizeof(status),
                 "%s  |  %d/%d done   survived %d   collapsed %d (%.0f%%)   |  %d workers",
                 state, completed, total, survived, collapsed, rate, workers);
    } else {
        snprintf(status, sizeof(status), "Configure a batch and press Run >>");
    }
    DrawText(status, 8, 40, 14, running ? (Color){0,230,120,255}
                                        : stopping ? (Color){255,180,0,255} : LIGHTGRAY);

    // progress bar
    if (d->has_batch && total > 0) {
        float fr = (float)completed / total;
        DrawRectangle(8, 24, plot_area_w - 8, 10, (Color){30,40,30,255});
        DrawRectangle(8, 24, (int)((plot_area_w - 8) * fr), 10, BAR);
    }

    // ----- control panel (right) -----
    DrawRectangle(panel_x, 0, panel_w, screen_h, (Color){25, 25, 30, 255});
    DrawLine(panel_x, 0, panel_x, screen_h, GRID);
    DrawText("Batch Setup", panel_x + 12, 12, 22, WHITE);

    float y = 50;
    float sw = panel_w - 200;  // slider width

    // Essentials
    slider_row(panel_x, y, sw, "Trajectories", TextFormat("%d", (int)d->f_num_traj), &d->f_num_traj, 10, 500); y += 30;
    slider_row(panel_x, y, sw, "Start size", TextFormat("%d", (int)d->f_initial), &d->f_initial, 1000, 20000); y += 30;
    slider_row(panel_x, y, sw, "Generations", TextFormat("%.1fM", d->f_target_gens/1e6f), &d->f_target_gens, 100000, 6000000); y += 32;
    GuiCheckBox((Rectangle){panel_x + 12, y, 20, 20}, "Elastic bounds", &d->elastic); y += 30;
    if (d->elastic) {
        slider_sweep_row(panel_x, y, sw, "Bound strength", TextFormat("%.2f", d->f_elasticity),
                   &d->f_elasticity, 0.0f, 1.0f, &d->sweep_enabled[5]); y += 30;
    }
    y += 4;

    bool busy = d->has_batch && (running || stopping);
    bool sweep_selected = d->sweep_enabled[0] || d->sweep_enabled[1] || d->sweep_enabled[2] || d->sweep_enabled[3] || d->sweep_enabled[4] || d->sweep_enabled[5];

    if (d->sweep_running) {
        float sweep_frac = (float)d->sweep_current_run / d->sweep_total_runs;
        DrawRectangle(panel_x + 12, y, (int)((panel_w - 24) * sweep_frac), 36, (Color){100, 150, 100, 255});
        if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36}, TextFormat("#133# Stop Sweep (%d/%d)", d->sweep_current_run, d->sweep_total_runs))) {
            d->sweep_running = false;
            batch_request_stop(&d->batch);
        }
    } else if (stopping) {
        GuiDisable();
        GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36}, "#133# Stopping...");
        GuiEnable();
    } else if (busy) {
        if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36}, "#133# Stop"))
            batch_request_stop(&d->batch);
    } else {
        if (sweep_selected) {
            if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36}, "#131# Start Sweep")) sweep_start(d);
        } else {
            if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36}, "#131# Run batch")) launch(d);
        }
    }
    y += 46;

    // Export the current histogram data to a timestamped CSV (cwd). Disabled until a
    // batch exists; shows the written filename briefly after a successful save.
    static char export_msg[280] = {0};
    static double export_at = 0.0;
    if (!d->has_batch) GuiDisable();
    if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 28}, "#02# Export data (CSV)")) {
        char path[256];
        export_histograms(d, NULL, path, sizeof(path));
        if (path[0]) snprintf(export_msg, sizeof(export_msg), "Saved %s", path);
        else         snprintf(export_msg, sizeof(export_msg), "Export failed");
        export_at = GetTime();
    }
    if (!d->has_batch) GuiEnable();
    y += 32;
    if (export_msg[0] && GetTime() - export_at < 6.0) {
        DrawText(export_msg, panel_x + 14, (int)y, 10, (Color){90, 200, 130, 255}); y += 14;
    }

    // Advanced (collapsible)
    if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 24},
                  d->show_advanced ? "#118# Advanced Options" : "#117# Advanced Options"))
        d->show_advanced = !d->show_advanced;
    y += 30;
    if (d->show_advanced) {
        slider_sweep_row(panel_x, y, sw, "Collapse <", TextFormat("%d", (int)d->f_collapse), &d->f_collapse, 0, 5000, &d->sweep_enabled[4]); y += 30;
        DrawText("Mutation", panel_x + 12, (int)y, 14, GRAY); y += 22;
        slider_sweep_row(panel_x, y, sw, "INDEL rate", TextFormat("%.2f", d->f_indel_rate), &d->f_indel_rate, 0.0f, 3.0f, &d->sweep_enabled[0]); y += 30;
        slider_sweep_row(panel_x, y, sw, "Mean size", TextFormat("%.1f", d->f_indel_size), &d->f_indel_size, 1.0f, 100.0f, &d->sweep_enabled[2]); y += 30;
        // Dup/del ratio: single dup:del SIZE ratio r (log slider centred at 1.0,
        // r in [0.001,1000]). Frequency is coupled as dup_bias = 1/(1+r) so the mean
        // array length stays constant (bigger events => proportionally rarer). The
        // four derived quantities (dup/del size and rate) are shown beneath.
        float ratio_e = log10f(d->f_size_ratio);
        const char *rfmt = d->f_size_ratio < 1.0f ? "%.3fx" : d->f_size_ratio < 10.0f ? "%.2fx" : "%.0fx";
        slider_sweep_row(panel_x, y, sw, "Dup/del ratio", TextFormat(rfmt, d->f_size_ratio), &ratio_e, -3.0f, 3.0f, &d->sweep_enabled[3]);
        d->f_size_ratio = powf(10.0f, ratio_e); y += 28;
        {
            float r = d->f_size_ratio, sr = sqrtf(r);
            float pdup = 1.0f / (1.0f + r);
            DrawText(TextFormat("dup ~%.1f u @ %.2f/gen", d->f_indel_size * sr, d->f_indel_rate * pdup),
                     panel_x + 14, (int)y, 10, (Color){90, 120, 90, 255}); y += 13;
            DrawText(TextFormat("del ~%.1f u @ %.2f/gen", d->f_indel_size / sr, d->f_indel_rate * (1.0f - pdup)),
                     panel_x + 14, (int)y, 10, (Color){90, 120, 90, 255}); y += 17;
        }
        slider_sweep_row(panel_x, y, sw, "SNP rate", TextFormat("%.2f", d->f_snp_rate), &d->f_snp_rate, 0.0f, 1.0f, &d->sweep_enabled[1]); y += 32;
        
        DrawText("Sweep", panel_x + 12, (int)y, 14, GRAY); y += 22;
        slider_row(panel_x, y, sw, "Sweep steps", TextFormat("%d", (int)d->f_sweep_steps), &d->f_sweep_steps, 2.0f, 50.0f); y += 30;
        slider_row(panel_x, y, sw, "Max run time (min)", TextFormat("%.1f", d->f_sweep_max_min), &d->f_sweep_max_min, 0.5f, 30.0f); y += 32;
        
        DrawText("Display", panel_x + 12, (int)y, 14, GRAY); y += 22;
        slider_row(panel_x, y, sw, "Display bars", TextFormat("%d", (int)d->f_nbins), &d->f_nbins, 16, 120); y += 26;
        DrawText("(adaptive; updates live)", panel_x + 14, (int)y, 10, (Color){90,120,90,255}); y += 18;
        slider_row(panel_x, y, sw, "Chebyshev order", TextFormat("%d", (int)d->f_cheb_order), &d->f_cheb_order, 1.0f, 9.0f); y += 26;
        GuiCheckBox((Rectangle){panel_x + 12, y, 20, 20}, "Autoscale X (fit data)", &d->autoscale_x); y += 26;
        GuiCheckBox((Rectangle){panel_x + 12, y, 20, 20}, "Show real data", &d->show_ref); y += 26;
        DrawText("Click \"Y:log/lin\" on a plot to toggle", panel_x + 14, (int)y, 10, GRAY); y += 22;
    }
    y += 6;

    // legend
    DrawText("overlays:", panel_x + 12, (int)y, 11, GRAY); y += 16;
    DrawRectangle(panel_x + 12, (int)y + 2, 12, 8, REAL);
    if (d->ref.loaded)
        DrawText(TextFormat("real (n=%d arrays)", d->ref.narrays), panel_x + 30, (int)y, 11, LIGHTGRAY);
    else
        DrawText("real data (no reference.dat)", panel_x + 30, (int)y, 11, (Color){150,120,60,255});
    y += 16;
    DrawRectangle(panel_x + 12, (int)y + 2, 12, 8, MED);  DrawText("live median of data", panel_x + 30, (int)y, 11, LIGHTGRAY); y += 16;

    y += 10;
    DrawText("Fit Parameters", panel_x + 12, (int)y, 14, GRAY); y += 22;
    for (int i = 0; i < 6; i++) {
        if (d->fit_text[i][0] != '\0') {
            DrawText(plots[i].title, panel_x + 14, (int)y, 11, LIGHTGRAY); y += 14;
            DrawText(d->fit_text[i], panel_x + 14, (int)y, 10, (Color){90, 120, 90, 255}); y += 16;
        }
    }

    // Sweep log intercept
    if (d->sweep_running && d->has_batch) {
        float elapsed = GetTime() - d->sweep_batch_start_time;
        if (complete || elapsed > d->f_sweep_max_min * 60.0) {
            sweep_log_and_advance(d, &su, &sh, &bs, &bg, &si, &cp);
        }
    }
}

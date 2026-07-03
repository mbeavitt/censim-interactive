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

static float g_gdens[DASH_MAXBINS];  // ghost per-bar PROBABILITY; reused each frame

// Aggregate the ghost histogram into `bars` display bars as a PROBABILITY (fraction
// of the ghost dataset) on the display axis, matching the sim's probability bars.
// Returns the peak so the sim and ghost can share one y-scale. Fills g_gdens.
static float ghost_fill_prob(const Histogram *ref, float lo_v, float hi_v, int log_x,
                             int bars, bool count_y) {
    if (!ref || ref->total <= 0 || bars < 1) return 0.0f;
    if (bars > DASH_MAXBINS) bars = DASH_MAXBINS;
    for (int j = 0; j < bars; j++) g_gdens[j] = 0.0f;
    for (int i = 0; i < ref->nbins; i++) {
        if (ref->counts[i] == 0) continue;
        float c = ref->log_scale ? sqrtf(edge_at(ref, i) * edge_at(ref, i + 1))
                                 : 0.5f * (edge_at(ref, i) + edge_at(ref, i + 1));
        float f = val_frac(c, lo_v, hi_v, log_x);
        if (f < 0.0f || f >= 1.0f) continue;
        int j = (int)(f * bars); if (j >= bars) j = bars - 1;
        g_gdens[j] += (float)ref->counts[i];
    }
    float peak = 0.0f;
    float L = log10f(lo_v > 0 ? lo_v : 1e-6f), H = log10f(hi_v > 0 ? hi_v : 1e-6f);
    for (int j = 0; j < bars; j++) {
        if (!count_y) {
            float vlo = log_x ? powf(10.0f, L + (H-L)*j/bars)     : lo_v + (hi_v-lo_v)*j/bars;
            float vhi = log_x ? powf(10.0f, L + (H-L)*(j+1)/bars) : lo_v + (hi_v-lo_v)*(j+1)/bars;
            float w = vhi - vlo; if (w <= 0.0f) w = 1.0f;
            g_gdens[j] /= w;             // density
        }
        g_gdens[j] /= (float)ref->total; // -> probability (fraction of dataset)
        if (g_gdens[j] > peak) peak = g_gdens[j];
    }
    return peak;
}

// Draw the ghost (faded bars + fit) from the pre-filled g_gdens probabilities, on
// the SHARED y-scale (dmax / lminv / lmaxv) so sim and ghost share one axis.
static void draw_ghost_overlay(const Histogram *ref, float px, float py, float pw, float ph,
                               float lo_v, float hi_v, int log_x, bool log_y, int bars,
                               int fit_type, int poly_order, bool do_bars, bool do_fit,
                               float bar_alpha, bool count_y,
                               float dmax, float lminv, float lmaxv,
                               char *out_fit_text) {
    if (out_fit_text) out_fit_text[0] = '\0';
    if (!ref || ref->total == 0 || bars < 1 || dmax <= 0.0f) return;
    if (bars > DASH_MAXBINS) bars = DASH_MAXBINS;

    #define GVAT(fr) (log_x ? powf(10.0f, log10f(lo_v) + (log10f(hi_v) - log10f(lo_v)) * (fr)) \
                            : lo_v + (hi_v - lo_v) * (fr))
    #define GH01(d) (log_y ? (log10f(d) - lminv) / (lmaxv - lminv) : (d) / dmax)

    // Faded bars.
    if (do_bars) {
        int a = (int)(bar_alpha * 255.0f + 0.5f); if (a < 0) a = 0; else if (a > 255) a = 255;
        Color fill = (Color){ REAL.r, REAL.g, REAL.b, (unsigned char)a };
        for (int j = 0; j < bars; j++) {
            if (g_gdens[j] <= 0.0f) continue;
            float h01 = GH01(g_gdens[j]); if (h01 < 0) h01 = 0; else if (h01 > 1) h01 = 1;
            float hh = ph * h01; if (log_y && hh < 1.0f) hh = 1.0f;
            int xL = (int)(px + pw * (float)j / bars);
            int xR = (int)(px + pw * (float)(j + 1) / bars);
            int w = xR - xL - 1; if (w < 1) w = 1;
            DrawRectangle(xL, (int)(py + ph - hh), w, (int)hh, fill);
        }
    }
    if (!do_fit) { return; }

    // ---- fit the ghost bars, same fit_type/order as the batch ----
    if (fit_type == 1 || fit_type == 3) {          // Chebyshev / power law (log-density)
        int n_coeffs = ((fit_type == 1) ? poly_order : 1) + 1;
        if (n_coeffs > 10) n_coeffs = 10;
        double M[100] = {0}, rhs[10] = {0}; int nfit = 0;
        double umin = log(lo_v > 0 ? lo_v : 1e-6), umax = log(hi_v > 0 ? hi_v : 1e-6);
        double fu_lo = 1e300, fu_hi = -1e300;
        for (int j = 0; j < bars; j++) {
            if (g_gdens[j] <= 0.0f) continue;
            float vlo = GVAT((float)j / bars), vhi = GVAT((float)(j + 1) / bars);
            double xc = log_x ? sqrt((double)vlo * vhi) : 0.5 * (vlo + vhi);
            if (xc <= 0) continue;
            double u = log(xc), y = log((double)g_gdens[j]);
            if (u < fu_lo) fu_lo = u; if (u > fu_hi) fu_hi = u;
            double z = (umax > umin) ? 2.0 * (u - umin) / (umax - umin) - 1.0 : 0.0;
            double T[10]; T[0] = 1; T[1] = z;
            for (int k = 2; k < n_coeffs; k++) T[k] = 2.0 * z * T[k-1] - T[k-2];
            for (int r = 0; r < n_coeffs; r++) {
                for (int cc = 0; cc < n_coeffs; cc++) M[r*n_coeffs+cc] += T[r]*T[cc];
                rhs[r] += T[r] * y;
            }
            nfit++;
        }
        double cf[10];
        if (nfit >= n_coeffs && solve_sys(M, rhs, cf, n_coeffs)) {
            if (out_fit_text) {   // same readout format as the sim fit
                double sse = 0.0;
                for (int j = 0; j < bars; j++) {
                    if (g_gdens[j] <= 0.0f) continue;
                    float vlo = GVAT((float)j / bars), vhi = GVAT((float)(j + 1) / bars);
                    double xc = log_x ? sqrt((double)vlo * vhi) : 0.5 * (vlo + vhi);
                    if (xc <= 0) continue;
                    double u = log(xc), y = log((double)g_gdens[j]);
                    double z = (umax > umin) ? 2.0 * (u - umin) / (umax - umin) - 1.0 : 0.0;
                    double T[10]; T[0] = 1; T[1] = z;
                    for (int k = 2; k < n_coeffs; k++) T[k] = 2.0 * z * T[k-1] - T[k-2];
                    double pred = 0; for (int k = 0; k < n_coeffs; k++) pred += cf[k]*T[k];
                    sse += (y - pred) * (y - pred);
                }
                int pos = 0;
                for (int i = 0; i < n_coeffs && pos < 120; i++)
                    pos += snprintf(out_fit_text + pos, 128 - pos, "c%d=%.2g ", i, cf[i]);
                snprintf(out_fit_text + pos, 128 - pos, "sse=%.2g", sse);
            }
            double pad = (nfit > 1) ? (fu_hi - fu_lo) / (double)(nfit - 1) : 0.0;
            double dlo = fu_lo - pad, dhi = fu_hi + pad;
            int hv = 0; float xprev = 0, yprev = 0;
            for (int k = 0; k <= 64; k++) {
                float fr = (float)k / 64.0f;
                float xv = GVAT(fr);
                if (xv <= 0) continue;
                double u = log(xv);
                if (u < dlo || u > dhi) { hv = 0; continue; }
                double z = (umax > umin) ? 2.0 * (u - umin) / (umax - umin) - 1.0 : 0.0;
                double T[10]; T[0] = 1; T[1] = z;
                for (int i = 2; i < n_coeffs; i++) T[i] = 2.0 * z * T[i-1] - T[i-2];
                double ld = 0; for (int i = 0; i < n_coeffs; i++) ld += cf[i]*T[i];
                double dd = exp(ld); if (!(dd > 0)) { hv = 0; continue; }
                float h01 = GH01((float)dd); if (h01 < 0) h01 = 0; else if (h01 > 1) h01 = 1;
                float X = px + fr * pw, Y = py + ph * (1.0f - h01);
                if (hv) DrawLineEx((Vector2){xprev, yprev}, (Vector2){X, Y}, 2.0f, REAL);
                xprev = X; yprev = Y; hv = 1;
            }
        }
    } else if (fit_type == 2 || fit_type == 4) {   // normal / gamma (from fine bins)
        double sx = 0, sx2 = 0; long tot = 0;
        for (int i = 0; i < ref->nbins; i++) {
            if (ref->counts[i] <= 0) continue;
            double x = ref->log_scale ? sqrt((double)edge_at(ref, i) * edge_at(ref, i + 1))
                                      : 0.5 * (edge_at(ref, i) + edge_at(ref, i + 1));
            sx += x * ref->counts[i]; sx2 += x * x * ref->counts[i]; tot += ref->counts[i];
        }
        if (tot > 1) {
            double mu = sx / tot, var = (sx2 - sx*sx/tot) / (tot - 1);
            double sd = var > 0 ? sqrt(var) : 0;
            if (sd > 0 && mu > 0) {
                double k_sh = (mu*mu)/(sd*sd), th = (sd*sd)/mu;   // gamma params
                #define GPDF(xv) (fit_type==2 \
                    ? exp(-0.5*(((xv)-mu)/sd)*(((xv)-mu)/sd)) / (sd*sqrt(2*M_PI)) \
                    : ((xv)>0 && k_sh>0 ? pow((xv),k_sh-1)*exp(-(xv)/th)/(pow(th,k_sh)*tgamma(k_sh)) : 0))
                // Probability: fraction-per-bar = pdf*width (count mode) or pdf (density),
                // on the same shared scale as the bars -- no peak re-fit needed.
                #define GWLOC(fr) (count_y ? (double)(GVAT(fminf((fr)+0.5f/bars,1.0f)) - GVAT(fmaxf((fr)-0.5f/bars,0.0f))) : 1.0)
                {
                    int hv = 0; float xp = 0, yp = 0;
                    for (int kk = 0; kk <= 64; kk++) {
                        float fr = (float)kk/64.0f; float xv = GVAT(fr);
                        double dd = GPDF(xv) * GWLOC(fr);
                        if (!(dd > 0)) { hv = 0; continue; }
                        float h01 = GH01((float)dd); if (h01<0) h01=0; else if (h01>1) h01=1;
                        float X = px + fr*pw, Y = py + ph*(1.0f - h01);
                        if (hv) DrawLineEx((Vector2){xp,yp},(Vector2){X,Y},2.0f, REAL);
                        xp = X; yp = Y; hv = 1;
                    }
                }
                if (out_fit_text) {   // same readout format as the sim fit
                    double sse = 0.0;
                    for (int j = 0; j < bars; j++) {
                        if (g_gdens[j] <= 0.0f) continue;
                        float fr = ((float)j + 0.5f) / bars;
                        double pred = GPDF(GVAT(fr)) * GWLOC(fr);
                        double diff = (double)g_gdens[j] - pred;
                        sse += diff * diff;
                    }
                    if (fit_type == 2) snprintf(out_fit_text, 128, "mu=%.2g, std=%.2g, sse=%.2g", mu, sd, sse);
                    else snprintf(out_fit_text, 128, "k=%.2g, theta=%.2g, sse=%.2g", k_sh, th, sse);
                }
                #undef GPDF
                #undef GWLOC
            }
        }
    }
    #undef GVAT
    #undef GH01
}

static void draw_hist(Rectangle b, const char *title, const HistSnap *s,
                      bool log_y, bool autoscale, const Histogram *ref,
                      int target_bars, Color accent, int fit_type, int poly_order, bool log_x,
                      bool count_y, bool ghost_bars, bool ghost_fit, float ghost_alpha,
                      char *out_fit_text, char *out_ghost_fit_text) {
    if (out_fit_text) out_fit_text[0] = '\0';
    if (out_ghost_fit_text) out_ghost_fit_text[0] = '\0';
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

    // A log-X display can't place a zero/negative lower bound (linear metrics like
    // similarity or HORs/kb start at 0), which would collapse the whole axis. Floor
    // lo_v to the first positive bin edge so log-X still works for those.
    if (log_x && lo_v <= 0.0f) {
        lo_v = bin_edge(s, b0 + 1);
        if (lo_v <= 0.0f) lo_v = (hi_v > 0.0f) ? hi_v * 1e-4f : 1e-6f;
    }

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
    // log-x span of the bars that actually constrain the fit; the overlay is
    // clamped to this so a high-order polynomial can't extrapolate (and blow up
    // via exp()) into the empty tail when the window is widened for the reference.
    double fit_u_lo = 1e300, fit_u_hi = -1e300;
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

    // Aggregate the window's fine bins into `bars` display bars BY VALUE on the
    // chosen display axis (log_x), so the x-scale can be toggled independently of
    // how the data was binned. Density = counts / linear bar width; the display
    // bars' value centres feed the fits below.
    static double bar_xc[DASH_MAXBINS];   // single-threaded UI: reused each frame
    static double bar_w[DASH_MAXBINS];    // value width of each display bar
    (void)avail;
    #define DVAT(fr) (log_x ? powf(10.0f, log10f(lo_v) + (log10f(hi_v) - log10f(lo_v)) * (fr)) \
                            : lo_v + (hi_v - lo_v) * (fr))
    for (int j = 0; j < bars; j++) dens[j] = 0.0f;
    for (int i = b0; i <= b1; i++) {
        if (s->counts[i] == 0) continue;
        float cen = s->log_scale ? sqrtf(bin_edge(s, i) * bin_edge(s, i + 1))
                                 : 0.5f * (bin_edge(s, i) + bin_edge(s, i + 1));
        float fr = val_frac(cen, lo_v, hi_v, log_x);
        if (fr < 0.0f || fr >= 1.0f) continue;
        int j = (int)(fr * bars); if (j >= bars) j = bars - 1;
        dens[j] += (float)s->counts[i];
    }
    for (int j = 0; j < bars; j++) {
        float vlo = DVAT((float)j / bars), vhi = DVAT((float)(j + 1) / bars);
        float w = vhi - vlo; if (w <= 0.0f) w = 1.0f;
        bar_w[j] = w;
        // count_y (per-array metrics): fraction of sims in this bin; otherwise a
        // probability density (fraction per unit x). Dividing by total puts the sim
        // and the ghost on ONE shared probability scale.
        if (!count_y) dens[j] /= w;
        if (s->total > 0) dens[j] /= (float)s->total;
        bar_xc[j] = log_x ? sqrt((double)vlo * vhi) : 0.5 * (vlo + vhi);
        float v = dens[j];
        if (v > 0.0f) {
            if (!seen || v > dmax_v) dmax_v = v;
            if (!seen || v < dmin_v) dmin_v = v;
            seen = 1;
        }
        if ((fit_type == 1 || fit_type == 3) && v > 0.0f && bar_xc[j] > 0.0) {
            double u = log(bar_xc[j]), y = log((double)v);
            if (u < fit_u_lo) fit_u_lo = u;
            if (u > fit_u_hi) fit_u_hi = u;
            double z = (u_max > u_min) ? 2.0 * (u - u_min) / (u_max - u_min) - 1.0 : 0.0;
            double T[10]; T[0] = 1.0; T[1] = z;
            for (int k = 2; k < n_coeffs; k++) T[k] = 2.0 * z * T[k-1] - T[k-2];
            for (int row = 0; row < n_coeffs; row++) {
                for (int col = 0; col < n_coeffs; col++) M_cheb[row*n_coeffs+col] += T[row]*T[col];
                rhs_cheb[row] += T[row] * y;
            }
            nfit++;
        }
    }
    // Fold the ghost's probability peak into the shared y-scale (also fills g_gdens
    // used when the ghost is drawn below), so sim + ghost sit on ONE axis.
    float ghost_peak = ghost_fill_prob(ref, lo_v, hi_v, log_x, bars, count_y);
    if (ghost_peak > dmax_v) dmax_v = ghost_peak;

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

    // y-axis: probability scale shared by sim + ghost. count_y shows % of the
    // dataset in the bin; otherwise a probability density (fraction per unit x).
    const char *ytop = count_y ? TextFormat("%.0f%%", dmax_v * 100.0f) : fmt(dmax_v);
    DrawText(ytop, (int)b.x + 3, (int)py - 4, 9, axc);
    DrawText(log_y ? (count_y ? TextFormat("%.1f%%", dmin_v*100.0f) : fmt(dmin_v)) : "0",
             (int)b.x + 3, (int)(py + ph - 8), 9, axc);
    DrawText(log_y ? "log" : "lin", (int)b.x + 3, (int)(py + ph/2), 9, axc);
    DrawText(count_y ? "%" : "pdf", (int)b.x + 3, (int)(py + ph/2 + 11), 8, shade(accent, 0.3f));

    // x-axis: 3 ticks across the displayed window (geometric mid for log)
    float xmid = log_x ? sqrtf(lo_v * hi_v) : 0.5f * (lo_v + hi_v);
    int ytick = (int)(py + ph + 2);
    DrawText(fmt(lo_v), (int)px, ytick, 9, axc);
    const char *m = fmt(xmid); DrawText(m, (int)(px + pw/2 - MeasureText(m, 9)/2), ytick, 9, axc);
    const char *x = fmt(hi_v); DrawText(x, (int)(px + pw - MeasureText(x, 9)), ytick, 9, axc);
    const char *xs = log_x ? "log x" : "lin x";
    DrawText(xs, (int)(px + pw - MeasureText(xs, 9)), (int)py - 2, 9, shade(accent, 0.3f));

    // ghost overlay (faded bars + fit), on the same display axis
    draw_ghost_overlay(ref, px, py, pw, ph, lo_v, hi_v, log_x, log_y, bars,
                       fit_type, poly_order, ghost_bars, ghost_fit, ghost_alpha, count_y,
                       dmax_v, lminv, lmaxv, out_ghost_fit_text);

    // Polynomial overlays (Chebyshev or Power law)
    if ((fit_type == 1 || fit_type == 3) && nfit >= n_coeffs) {
        double cf[10];
        if (solve_sys(M_cheb, rhs_cheb, cf, n_coeffs)) {
            double sse = 0.0;
            for (int j = 0; j < bars; j++) {
                if (dens[j] <= 0.0f) continue;
                double xc = bar_xc[j];
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
            // Allow one extra bar-spacing beyond the outermost fitted bars so the
            // curve runs down to the axis instead of stopping half a bar short.
            double pad = (nfit > 1) ? (fit_u_hi - fit_u_lo) / (double)(nfit - 1) : 0.0;
            double draw_lo = fit_u_lo - pad, draw_hi = fit_u_hi + pad;
            int hv = 0; float xprev = 0, yprev = 0;
            for (int k = 0; k <= 64; k++) {
                float f = (float)k / 64.0f;
                float xv = DVAT(f);
                if (xv <= 0.0) continue;
                double u = log(xv);
                // Don't extrapolate far past the fitted data: break the line outside
                // the populated range (+1 bar) so it can't ping back up in the tail.
                if (u < draw_lo || u > draw_hi) { hv = 0; continue; }
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
        double theta = (stddev * stddev) / mu;
        double k_shape = (mu * mu) / (stddev * stddev);
        for (int k = 0; k <= 64; k++) {
            float f = (float)k / 64.0f;
            float xv = DVAT(f);
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
            // probability: pdf (density) or pdf*width (count mode: fraction per bar)
            double wloc = count_y ? (DVAT(fminf(f + 0.5f/bars, 1.0f)) - DVAT(fmaxf(f - 0.5f/bars, 0.0f))) : 1.0;
            double d = pdf_x * wloc;
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
            double xc = bar_xc[j];
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
            double d_pred = pdf_x * (count_y ? bar_w[j] : 1.0);
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
        float f = val_frac(med, lo_v, hi_v, log_x);
        if (f >= 0 && f <= 1) {
            float xx = px + f * pw;
            for (int yy = (int)py; yy < (int)(py + ph); yy += 7)
                DrawLineEx((Vector2){xx, (float)yy}, (Vector2){xx, (float)yy + 4}, 2.5f, MED);
            char mb[32]; snprintf(mb, sizeof(mb), "med %s", fmt(med));
            DrawText(mb, (int)xx + 2, (int)(py + ph - 11), 9, MED);
        }
    }
    #undef DVAT
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
        d->ghost_fit_text[i][0] = '\0';
    }
    d->autoscale_x = false;
    // Per-plot log-Y defaults (left->right, top->bottom): lin lin log / log lin log
    d->plot_log_y[0] = false;  // unique/kb   lin
    d->plot_log_y[1] = false;  // HORs/kb     lin
    d->plot_log_y[2] = true;   // block size  log
    d->plot_log_y[3] = true;   // block gap   log
    d->plot_log_y[4] = false;  // similarity  lin
    d->plot_log_y[5] = true;   // composite   log
    // Per-plot log-X defaults: match each metric's natural (binned) scale.
    d->plot_log_x[0] = false;  // unique/kb   lin
    d->plot_log_x[1] = true;   // HORs/kb     log (spans orders of magnitude)
    d->plot_log_x[2] = true;   // block size  log
    d->plot_log_x[3] = true;   // block gap   log
    d->plot_log_x[4] = false;  // similarity  lin
    d->plot_log_x[5] = true;   // composite   log

    // Default ghost: the real Arabidopsis histogram (66 genomes). Try $CENSIM_REFERENCE,
    // ./data/... , then next to / above the executable, then the legacy reference.dat.
    reference_init(&d->ref);
    d->show_ref = true;
    d->ghost_bars = true;   // faded histogram bars
    d->ghost_fit = true;    // fitted curve, same fit as the batch
    d->ghost_alpha = 0.28f; // ghost bar opacity

    const char *env = getenv("CENSIM_REFERENCE");
    const char *appdir = GetApplicationDirectory();
    char p1[1024], p2[1024];
    snprintf(p1, sizeof(p1), "%sdata/real_arabidopsis_hist.csv", appdir);
    snprintf(p2, sizeof(p2), "%s../../data/real_arabidopsis_hist.csv", appdir);
    bool got = (env && reference_load_run(&d->ref, env))
             || reference_load_run(&d->ref, "data/real_arabidopsis_hist.csv")
             || reference_load_run(&d->ref, p1)
             || reference_load_run(&d->ref, p2);
    if (got) {
        snprintf(d->ghost_name, sizeof(d->ghost_name), "real Arabidopsis (66 genomes)");
    } else if (!reference_load(&d->ref, "reference.dat")) {   // legacy fallback
        char path[1024];
        snprintf(path, sizeof(path), "%sreference.dat", appdir);
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

static int import_histograms(Dashboard *d, const char *path) {
    if (d->has_batch) {
        batch_request_stop(&d->batch);
        batch_join(&d->batch);
        batch_free(&d->batch);
        d->has_batch = false;
    }
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    
    memset(&d->batch, 0, sizeof(Batch));
    pthread_mutex_init(&d->batch.lock, NULL);
    
    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') {
            if (strncmp(line, "# trajectories=", 15) == 0) {
                int traj, surv, coll;
                if (sscanf(line, "# trajectories=%d survived=%d collapsed=%d", &traj, &surv, &coll) == 3) {
                    d->batch.cfg.num_trajectories = traj;
                    d->batch.survived = surv;
                    d->batch.collapsed_count = coll;
                    d->batch.completed = traj;
                }
            } else {
                char metric[64];
                float min_v, max_v;
                int nbins, log_scale;
                long total, underflow, overflow;
                if (sscanf(line, "# %63[^:]: min=%f max=%f nbins=%d log_scale=%d total=%ld underflow=%ld overflow=%ld",
                           metric, &min_v, &max_v, &nbins, &log_scale, &total, &underflow, &overflow) == 8) {
                    Histogram *h = NULL;
                    if (strcmp(metric, "unique_per_kb") == 0) h = &d->batch.h_unique_per_kb;
                    else if (strcmp(metric, "hors_per_kb") == 0) h = &d->batch.h_hors_per_kb;
                    else if (strcmp(metric, "block_size") == 0) h = &d->batch.h_block_size;
                    else if (strcmp(metric, "block_gap") == 0) h = &d->batch.h_block_gap;
                    else if (strcmp(metric, "similarity") == 0) h = &d->batch.h_similarity;
                    else if (strcmp(metric, "diversity") == 0) h = &d->batch.h_diversity;
                    else if (strcmp(metric, "composite") == 0) h = &d->batch.h_composite;
                    else if (strcmp(metric, "collapse_gen") == 0) h = &d->batch.h_collapse_gen;
                    
                    if (h) {
                        hist_init(h, min_v, max_v, nbins, log_scale);
                        h->total = total;
                        h->underflow = underflow;
                        h->overflow = overflow;
                    }
                }
            }
        } else if (strncmp(line, "metric", 6) != 0) {
            char metric[64];
            int bin;
            long count;
            char *comma1 = strchr(line, ',');
            if (comma1) {
                int len = comma1 - line;
                if (len > 63) len = 63;
                strncpy(metric, line, len);
                metric[len] = 0;
                
                char *comma5 = strrchr(line, ',');
                if (comma5) {
                    count = atol(comma5 + 1);
                    bin = atoi(comma1 + 1);
                    
                    Histogram *h = NULL;
                    if (strcmp(metric, "unique_per_kb") == 0) h = &d->batch.h_unique_per_kb;
                    else if (strcmp(metric, "hors_per_kb") == 0) h = &d->batch.h_hors_per_kb;
                    else if (strcmp(metric, "block_size") == 0) h = &d->batch.h_block_size;
                    else if (strcmp(metric, "block_gap") == 0) h = &d->batch.h_block_gap;
                    else if (strcmp(metric, "similarity") == 0) h = &d->batch.h_similarity;
                    else if (strcmp(metric, "diversity") == 0) h = &d->batch.h_diversity;
                    else if (strcmp(metric, "composite") == 0) h = &d->batch.h_composite;
                    else if (strcmp(metric, "collapse_gen") == 0) h = &d->batch.h_collapse_gen;
                    
                    if (h && bin >= 0 && bin < h->nbins) {
                        h->counts[bin] = count;
                    }
                }
            }
        }
    }
    fclose(f);
    d->has_batch = true;
    d->started = true;
    d->batch.running = false;
    d->batch.stop_requested = false;
    d->batch.num_workers = 0;
    d->batch.workers_running = 0;
    return 1;
}

// A labelled slider row; returns the (possibly updated) value via the pointer.
// UI text scaling, shared with main.c's Options dialog. df() scales pixel
// offsets smoothly. df_font() scales font sizes but snaps to an integer multiple
// of the 10px default-font base -- the bitmap font is only crisp at 10/20/30...,
// so fractional sizes (e.g. 15px) point-upscale and look rough. It never shrinks
// below the authored size.
extern float g_ui_scale;
static int  df(int v) { int s = (int)(v * g_ui_scale + 0.5f); return s < 1 ? 1 : s; }
static int  df_font(int base) {
    int want = (int)(base * g_ui_scale + 0.5f);
    if (want < base) want = base;
    int snapped = ((want + 5) / 10) * 10;      // nearest multiple of 10
    if (snapped < base) snapped = base;        // don't shrink below the design size
    return snapped;
}
static void DrawTextS(const char *t, int x, int y, int sz, Color c) { DrawText(t, x, y, df_font(sz), c); }

// Lay out a slider so its right edge stays fixed as the (scaled) label grows:
// push the slider start right by df(130) and shrink its width to match.
static void slider_geom(float panel_x, float w, float *sx, float *sw) {
    *sx = panel_x + df(130);
    *sw = w + 130 - df(130);
    if (*sw < 40) *sw = 40;
}

static void slider_row(float panel_x, float y, float w, const char *label,
                       const char *valtext, float *v, float lo, float hi) {
    DrawTextS(label, (int)panel_x + 12, (int)y + 2, 14, LIGHTGRAY);
    float sx, sw; slider_geom(panel_x, w, &sx, &sw);
    GuiSlider((Rectangle){sx, y, sw, 18}, NULL, valtext, v, lo, hi);
}

static void slider_sweep_row(float panel_x, float y, float w, const char *label,
                             const char *valtext, float *v, float lo, float hi, bool *sweep) {
    if (sweep) {
        GuiCheckBox((Rectangle){panel_x + 6, y + 2, 14, 14}, "", sweep);
        DrawTextS(label, (int)panel_x + 26, (int)y + 2, 14, LIGHTGRAY);
    } else {
        DrawTextS(label, (int)panel_x + 12, (int)y + 2, 14, LIGHTGRAY);
    }
    float sx, sw; slider_geom(panel_x, w, &sx, &sw);
    GuiSlider((Rectangle){sx, y, sw, 18}, NULL, valtext, v, lo, hi);
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

// Load a sweep run's input parameters (1-based) from summary.csv so the
// browse-mode title reflects that run's actual inputs. Columns:
//   run,indel_rate,snp_rate,indel_size,size_ratio,collapse,elasticity,...
static void load_run_params(Dashboard *d, int run) {
    char path[512];
    snprintf(path, sizeof(path), "%s/summary.csv", d->sweep_dir);
    FILE *f = fopen(path, "r");
    if (!f) return;
    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        int r; float ir, sr, is, ratio, coll, elast;
        if (sscanf(line, "%d,%f,%f,%f,%f,%f,%f", &r, &ir, &sr, &is, &ratio, &coll, &elast) == 7
            && r == run) {
            d->f_indel_rate = ir;
            d->f_snp_rate   = sr;
            d->f_indel_size = is;
            d->f_size_ratio = ratio;
            d->f_collapse   = coll;
            d->f_elasticity = elast;
            d->elastic      = (elast > 0.0f);
            break;
        }
    }
    fclose(f);
}

// Open a native chooser (zenity on Linux, osascript on macOS) for a file or a
// directory, returning the selected path. Returns 1 on success / 0 if cancelled.
static int pick_path(const char *title, int directory, char *out, size_t n) {
    char cmd[512];
#ifdef __APPLE__
    snprintf(cmd, sizeof(cmd),
        "osascript -e 'POSIX path of (choose %s with prompt \"%s\")' 2>/dev/null",
        directory ? "folder" : "file", title);
#else
    snprintf(cmd, sizeof(cmd),
        "zenity --file-selection %s --title=\"%s\" 2>/dev/null",
        directory ? "--directory" : "", title);
#endif
    FILE *p = popen(cmd, "r");
    if (!p) return 0;
    int ok = 0;
    if (fgets(out, (int)n, p)) {
        out[strcspn(out, "\n")] = 0;
        if (out[0]) ok = 1;
    }
    pclose(p);
    return ok;
}

// Set the ghost overlay from a run/sweep histogram CSV and label it by filename.
static void set_ghost_from_file(Dashboard *d, const char *path) {
    if (reference_load_run(&d->ref, path)) {
        const char *base = strrchr(path, '/');
        snprintf(d->ghost_name, sizeof(d->ghost_name), "%s", base ? base + 1 : path);
        d->show_ref = true;
    }
}

// Count run_*.csv files in a sweep directory (for the browse scrubber).
static int count_runs(const char *dir) {
    char cmd[600];
    snprintf(cmd, sizeof(cmd), "ls \"%s\"/histograms/run_*.csv 2>/dev/null | wc -l", dir);
    FILE *p = popen(cmd, "r");
    int n = 0;
    if (p) { if (fscanf(p, "%d", &n) != 1) n = 0; pclose(p); }
    return n;
}

// Load browse run `idx` (0-based, clamped to the sweep), syncing slider + params.
static void browse_load(Dashboard *d, int idx) {
    if (idx < 0) idx = 0;
    if (d->browse_run_count > 0 && idx >= d->browse_run_count) idx = d->browse_run_count - 1;
    char path[512];
    snprintf(path, sizeof(path), "%s/histograms/run_%04d.csv", d->sweep_dir, idx + 1);
    if (import_histograms(d, path)) {
        d->browse_run_idx = idx;
        d->f_browse_run   = (float)(idx + 1);
        load_run_params(d, idx + 1);
    }
}

// Enter browse mode on a sweep directory (expects <dir>/histograms/run_XXXX.csv).
static void open_sweep_dir(Dashboard *d, const char *dir) {
    snprintf(d->sweep_dir, sizeof(d->sweep_dir), "%s", dir);
    d->browse_run_count = count_runs(dir);
    d->sweep_browsing = true;
    d->browse_run_idx = 0;
    d->f_browse_run = 1.0f;
    char path[512];
    snprintf(path, sizeof(path), "%s/histograms/run_0001.csv", d->sweep_dir);
    if (import_histograms(d, path)) load_run_params(d, 1);
    else d->sweep_browsing = false;   // not a sweep dir
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
                  ref, target_bars, accent, fit_type, (int)d->f_cheb_order, d->plot_log_x[i],
                  i < 2,   // unique/kb + HORs/kb are per-array: show frequency (sim count)
                  d->show_ref && d->ghost_bars, d->show_ref && d->ghost_fit, d->ghost_alpha,
                  d->fit_text[i], d->ghost_fit_text[i]);
        // per-plot clickable Y- and X-scale toggles (top-right of the cell)
        Rectangle tg = { cell[i].x + cell[i].width - 52, cell[i].y + 3, 46, 15 };
        bool hov = CheckCollisionPointRec(mouse, tg);
        DrawRectangleRec(tg, hov ? shade(accent, 0.3f) : (Color){25,35,25,255});
        DrawRectangleLinesEx(tg, 1, shade(accent, 0.55f));
        DrawText(d->plot_log_y[i] ? "Y:log" : "Y:lin", (int)tg.x + 5, (int)tg.y + 3, 9, accent);
        if (hov && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) d->plot_log_y[i] = !d->plot_log_y[i];

        Rectangle tx = { cell[i].x + cell[i].width - 102, cell[i].y + 3, 46, 15 };
        bool hovx = CheckCollisionPointRec(mouse, tx);
        DrawRectangleRec(tx, hovx ? shade(accent, 0.3f) : (Color){25,35,25,255});
        DrawRectangleLinesEx(tx, 1, shade(accent, 0.55f));
        DrawText(d->plot_log_x[i] ? "X:log" : "X:lin", (int)tx.x + 5, (int)tx.y + 3, 9, accent);
        if (hovx && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) d->plot_log_x[i] = !d->plot_log_x[i];
    }

    // ----- status bar -----
    DrawRectangle(0, top - 2, plot_area_w + 4, 2, GRID);
    
    // ----- title: current mode + input parameters (always visible) -----
    const char *mode = d->sweep_running  ? "SWEEP"
                     : d->sweep_browsing ? "BROWSING"
                     : d->has_batch      ? "LIVE"
                     :                     "READY";
    char title[360];
    snprintf(title, sizeof(title),
             "%s%s   |   INDEL rate %.2f  size %.1f   SNP %.2f   dup:del %.2fx   "
             "gens %.0f   init %.0f   collapse %.0f   elastic %s",
             mode,
             d->sweep_browsing ? TextFormat(" run %d", d->browse_run_idx + 1) : "",
             d->f_indel_rate, d->f_indel_size, d->f_snp_rate, d->f_size_ratio,
             d->f_target_gens, d->f_initial, d->f_collapse,
             d->elastic ? TextFormat("%.2f", d->f_elasticity) : "off");
    int tw = MeasureText(title, 14);
    int tx = (plot_area_w - tw) / 2; if (tx < 8) tx = 8;
    DrawText(title, tx, 40, 14, (Color){150, 200, 255, 255});   // centred, under the progress bar

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
    DrawText(status, 8, 8, 14, running ? (Color){0,230,120,255}
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

    // The panel content (below the fixed title) often overflows -- Advanced Options
    // and the Fit Parameters readout push well past the bottom edge. Let the wheel
    // scroll it, like the single view. Content height is measured each frame (the
    // final y below) and used to clamp the offset on the following frame.
    const int panel_head = 44;  // fixed title band; content scrolls beneath it
    if (mouse.x >= panel_x && mouse.x <= screen_w) {
        float wheel = GetMouseWheelMove();
        if (wheel != 0.0f) d->panel_scroll += wheel * 30.0f;
    }
    float max_scroll = -(d->panel_content_h - screen_h + 12);
    if (max_scroll > 0.0f) max_scroll = 0.0f;
    if (d->panel_scroll > 0.0f)         d->panel_scroll = 0.0f;
    if (d->panel_scroll < max_scroll)   d->panel_scroll = max_scroll;

    DrawText("Batch Setup", panel_x + 12, 12, 22, WHITE);

    BeginScissorMode(panel_x, panel_head, panel_w, screen_h - panel_head);

    float y = 50 + d->panel_scroll;
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
    } else if (d->sweep_browsing) {
        int total = d->browse_run_count > 0 ? d->browse_run_count : d->browse_run_idx + 1;
        if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36},
                      TextFormat("#131# Browse Mode (Run %d/%d)", d->browse_run_idx + 1, total))) {
            d->sweep_browsing = false;
        }
        // Scrub slider: drag across all runs; loads on each integer change.
        DrawText("Run", (int)panel_x + 12, (int)y + 44, 12, LIGHTGRAY);
        GuiSlider((Rectangle){panel_x + 46, y + 42, panel_w - 58, 18}, NULL,
                  TextFormat("%d", d->browse_run_idx + 1), &d->f_browse_run,
                  1.0f, (float)total);
        int want = (int)(d->f_browse_run + 0.5f) - 1;
        if (want != d->browse_run_idx) browse_load(d, want);

        if (GuiButton((Rectangle){panel_x + 12, y + 68, (panel_w - 28)/2, 28}, "#118# Prev"))
            browse_load(d, d->browse_run_idx - 1);
        if (GuiButton((Rectangle){panel_x + 16 + (panel_w - 28)/2, y + 68, (panel_w - 28)/2, 28}, "#119# Next"))
            browse_load(d, d->browse_run_idx + 1);

        // Pin the run you're looking at as the ghost, then scrub to another run to
        // compare it against -- past-run-vs-past-run overlay.
        if (GuiButton((Rectangle){panel_x + 12, y + 100, panel_w - 24, 24}, "#77# Pin this run as ghost")) {
            char path[512];
            snprintf(path, sizeof(path), "%s/histograms/run_%04d.csv", d->sweep_dir, d->browse_run_idx + 1);
            set_ghost_from_file(d, path);
        }
    } else {
        if (sweep_selected) {
            if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36}, "#131# Start Sweep")) sweep_start(d);
        } else {
            if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36}, "#131# Run batch")) launch(d);
        }
    }
    y += 46;
    if (d->sweep_browsing) y += 104;

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
    
    // Browser: latest sweep (quick) or pick any sweep directory.
    if (GuiButton((Rectangle){panel_x + 12, y, (panel_w - 28)/2, 28}, "#05# Latest")) {
        FILE *p = popen("ls -td sweep_* 2>/dev/null | head -n 1", "r");
        if (p) {
            char latest[256];
            if (fgets(latest, sizeof(latest), p)) {
                latest[strcspn(latest, "\n")] = 0;
                if (latest[0]) open_sweep_dir(d, latest);
            }
            pclose(p);
        }
    }
    if (GuiButton((Rectangle){panel_x + 16 + (panel_w - 28)/2, y, (panel_w - 28)/2, 28}, "#01# Open sweep...")) {
        char dir[256];
        if (pick_path("Choose a sweep directory", 1, dir, sizeof(dir))) open_sweep_dir(d, dir);
    }
    y += 32;

    // Ghost overlay: pick any run/sweep CSV (a past run, another sweep, real data)
    // to overlay on the current plots for direct comparison.
    if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 28}, "#146# Pick ghost overlay...")) {
        char file[512];
        if (pick_path("Choose a histogram CSV to overlay as ghost", 0, file, sizeof(file)))
            set_ghost_from_file(d, file);
    }
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
        GuiCheckBox((Rectangle){panel_x + 12, y, 20, 20}, "Show ghost overlay", &d->show_ref); y += 26;
        if (d->show_ref) {
            GuiCheckBox((Rectangle){panel_x + 28, y, 18, 18}, "ghost bars", &d->ghost_bars);
            GuiCheckBox((Rectangle){panel_x + 28 + (panel_w-40)/2, y, 18, 18}, "ghost fit", &d->ghost_fit);
            y += 24;
            if (d->ghost_bars) {
                slider_row(panel_x, y, sw, "Ghost opacity",
                           TextFormat("%.0f%%", d->ghost_alpha * 100.0f),
                           &d->ghost_alpha, 0.05f, 1.0f); y += 26;
            }
        }
        DrawText("Click \"Y:log/lin\" on a plot to toggle", panel_x + 14, (int)y, 10, GRAY); y += 22;
    }
    y += 6;

    // legend
    DrawTextS("overlays:", panel_x + 12, (int)y, 11, GRAY); y += df(16);
    DrawRectangle(panel_x + 12, (int)y + 2, 12, 8, REAL);
    if (d->ref.loaded) {
        const char *gname = d->ghost_name[0] ? d->ghost_name : "reference.dat";
        DrawTextS(TextFormat("ghost: %s", gname), panel_x + 30, (int)y, 11, LIGHTGRAY);
    } else {
        DrawTextS("ghost: none (Pick ghost overlay...)", panel_x + 30, (int)y, 11, (Color){150,120,60,255});
    }
    y += df(16);
    DrawRectangle(panel_x + 12, (int)y + 2, 12, 8, MED);  DrawTextS("live median of data", panel_x + 30, (int)y, 11, LIGHTGRAY); y += df(16);

    y += 10;
    DrawTextS("Sim Fit Parameters", panel_x + 12, (int)y, 14, GRAY); y += df(22);
    for (int i = 0; i < 6; i++) {
        if (d->fit_text[i][0] != '\0') {
            DrawTextS(plots[i].title, panel_x + 14, (int)y, 11, LIGHTGRAY); y += df(14);
            DrawTextS(d->fit_text[i], panel_x + 14, (int)y, 10, (Color){90, 120, 90, 255}); y += df(16);
        }
    }

    // Ghost fit parameters (only when a ghost fit is actually being drawn).
    // Must be drawn inside the scissor region and before panel_content_h is
    // recorded, otherwise it's excluded from the scroll clamp and never reachable.
    bool any_ghost_fit = false;
    for (int i = 0; i < 6; i++) if (d->ghost_fit_text[i][0] != '\0') { any_ghost_fit = true; break; }
    if (any_ghost_fit) {
        y += 10;
        DrawTextS("Ghost Fit Parameters", panel_x + 12, (int)y, 14, GRAY); y += df(22);
        for (int i = 0; i < 6; i++) {
            if (d->ghost_fit_text[i][0] != '\0') {
                DrawTextS(plots[i].title, panel_x + 14, (int)y, 11, LIGHTGRAY); y += df(14);
                DrawTextS(d->ghost_fit_text[i], panel_x + 14, (int)y, 10, (Color){130, 90, 90, 255}); y += df(16);
            }
        }
    }

    // Record the content bottom (in unscrolled coords) for next frame's clamp.
    d->panel_content_h = y - d->panel_scroll;
    EndScissorMode();

    // Sweep log intercept
    if (d->sweep_running && d->has_batch) {
        float elapsed = GetTime() - d->sweep_batch_start_time;
        if (complete || elapsed > d->f_sweep_max_min * 60.0) {
            sweep_log_and_advance(d, &su, &sh, &bs, &bg, &si, &cp);
        }
    }
}

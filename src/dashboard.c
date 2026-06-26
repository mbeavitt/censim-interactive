#include "dashboard.h"
#include "config.h"
#include "hist.h"
#include <raylib.h>
#include "raygui.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

// ---------------------------------------------------------------------------
// Reference ("ghost") values overlaid on the plots.
// TODO: use real data — these are STAND-INS estimated from the R1 report summary
// stats, not the empirical binned distributions. See repo TODO.md.
// ---------------------------------------------------------------------------
#define REF_REAL_UNIQUE_PER_KB     2.3f
#define REF_UNIFORM_UNIQUE_PER_KB  1.5f
#define REF_REAL_HORS_PER_KB       261.4f
#define REF_UNIFORM_HORS_PER_KB    124.5f
#define REF_REAL_BLOCK_SIZE_MED    3.5f
#define REF_REAL_SIMILARITY_MED    0.16f

#define DASH_MAXBINS 256

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

typedef struct {
    float value;
    Color color;
    const char *label;
} RefMark;

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

static const Color BG    = (Color){15, 20, 15, 255};
static const Color GRID  = (Color){0, 90, 40, 255};
static const Color BAR   = (Color){0, 200, 110, 255};
static const Color REAL  = (Color){255, 70, 70, 255};   // A. thaliana ghost
static const Color UNIF  = (Color){150, 150, 160, 255}; // uniform-model ghost (unused for now)
static const Color MED   = (Color){255, 230, 80, 255};  // live median of the plotted data

static void draw_hist(Rectangle b, const char *title, const HistSnap *s,
                      RefMark *refs, int nrefs, bool log_y, bool autoscale) {
    DrawRectangleRec(b, BG);
    DrawRectangleLinesEx(b, 1, GRID);
    DrawText(title, (int)b.x + 6, (int)b.y + 4, 11, BAR);
    char nbuf[40];
    snprintf(nbuf, sizeof(nbuf), "n=%ld", s ? s->total : 0);
    DrawText(nbuf, (int)b.x + 6, (int)b.y + 17, 9, GRID);

    // Warn if any data fell outside the binned range (clipped to under/overflow).
    if (s && (s->underflow > 0 || s->overflow > 0)) {
        char ob[48];
        snprintf(ob, sizeof(ob), "clipped <%ld >%ld", s->underflow, s->overflow);
        DrawText(ob, (int)(b.x + b.width) - MeasureText(ob, 9) - 6, (int)b.y + 17, 9, (Color){255,140,0,255});
    }

    // Inner plot area: leave a left margin for y numbers and a bottom strip for x.
    float px = b.x + 40, py = b.y + 30, pw = b.width - 48, ph = b.height - 46;

    if (!s || s->maxcount == 0) {
        DrawText("awaiting data", (int)(px + pw/2 - 40), (int)(py + ph/2), 10, GRID);
        return;
    }

    // displayed x-window (autoscaled to visible bins at this Y scale) and value range
    int b0, b1; window_bins(s, autoscale, log_y, &b0, &b1);
    // Always keep the natural-data reference lines in view, so the gap between the
    // simulated distribution and the real value stays visible even when autoscaled.
    for (int r = 0; r < nrefs; r++) {
        float rf = val_frac(refs[r].value, s->min, s->max, s->log_scale);
        if (rf < 0.0f) continue;
        int rb = (int)(rf * s->nbins);
        if (rb < 0) rb = 0; else if (rb >= s->nbins) rb = s->nbins - 1;
        if (rb < b0) b0 = rb;
        if (rb > b1) b1 = rb;
    }
    int nb = b1 - b0 + 1;
    float lo_v = bin_edge(s, b0), hi_v = bin_edge(s, b1 + 1);

    // y scales to the tallest *visible* bar
    long dmax = 0;
    for (int i = b0; i <= b1; i++) if (s->counts[i] > dmax) dmax = s->counts[i];
    if (dmax == 0) dmax = 1;

    // bars
    float denom = log_y ? log10f((float)dmax + 1.0f) : (float)dmax;
    float bw = pw / nb;
    for (int i = b0; i <= b1; i++) {
        if (s->counts[i] == 0) continue;
        float num = log_y ? log10f((float)s->counts[i] + 1.0f) : (float)s->counts[i];
        float hh = ph * num / denom;
        DrawRectangle((int)(px + (i - b0) * bw), (int)(py + ph - hh),
                      (int)(bw > 1 ? bw - 1 : 1), (int)hh, BAR);
    }

    // y-axis: tick at top (max visible count) and bottom (0); midline label
    DrawText(fmt((float)dmax), (int)b.x + 3, (int)py - 4, 9, GRID);
    DrawText("0", (int)b.x + 3, (int)(py + ph - 8), 9, GRID);
    DrawText(log_y ? "log" : "lin", (int)b.x + 3, (int)(py + ph/2), 9, GRID);

    // x-axis: 3 ticks across the displayed window (geometric mid for log)
    float xmid = s->log_scale ? sqrtf(lo_v * hi_v) : 0.5f * (lo_v + hi_v);
    int ytick = (int)(py + ph + 2);
    DrawText(fmt(lo_v), (int)px, ytick, 9, GRID);
    const char *m = fmt(xmid); DrawText(m, (int)(px + pw/2 - MeasureText(m, 9)/2), ytick, 9, GRID);
    const char *x = fmt(hi_v); DrawText(x, (int)(px + pw - MeasureText(x, 9)), ytick, 9, GRID);
    const char *xs = s->log_scale ? "log x" : "lin x";
    DrawText(xs, (int)(px + pw - MeasureText(xs, 9)), (int)py - 2, 9, (Color){0,70,30,255});

    // reference markers (vertical lines), positioned within the window
    for (int r = 0; r < nrefs; r++) {
        float f = val_frac(refs[r].value, lo_v, hi_v, s->log_scale);
        if (f < 0 || f > 1) continue;
        float xx = px + f * pw;
        DrawLine((int)xx, (int)py, (int)xx, (int)(py + ph), refs[r].color);
        if (refs[r].label) DrawText(refs[r].label, (int)xx + 2, (int)py + 1, 9, refs[r].color);
    }

    // live median of the plotted data (dashed bright marker)
    float med = snap_median(s);
    if (!isnan(med)) {
        float f = val_frac(med, lo_v, hi_v, s->log_scale);
        if (f >= 0 && f <= 1) {
            float xx = px + f * pw;
            for (int yy = (int)py; yy < (int)(py + ph); yy += 6)
                DrawLine((int)xx, yy, (int)xx, yy + 3, MED);
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
    d->unbounded     = true;
    d->f_nbins       = 50.0f;
    d->autoscale_x   = true;
    d->show_advanced = false;
    d->has_batch     = false;
    d->started       = false;
    // Per-plot log-Y defaults (left->right, top->bottom): lin lin log / lin log lin
    d->plot_log_y[0] = false;  // unique/kb   lin
    d->plot_log_y[1] = false;  // HORs/kb     lin
    d->plot_log_y[2] = true;   // block size  log
    d->plot_log_y[3] = false;  // block gap   lin
    d->plot_log_y[4] = true;   // similarity  log
    d->plot_log_y[5] = false;  // composite   lin
}

void dashboard_free(Dashboard *d) {
    if (d->has_batch) { batch_free(&d->batch); d->has_batch = false; }
}

static void launch(Dashboard *d) {
    if (d->has_batch) { batch_free(&d->batch); d->has_batch = false; }

    SimParams p;
    batch_default_params(&p, d->unbounded);
    p.indel_rate        = d->f_indel_rate;
    p.snp_rate          = d->f_snp_rate;
    p.indel_size_lambda = d->f_indel_size;
    p.collapse_threshold = (int)d->f_collapse;
    p.target_size       = (int)d->f_initial;  // elastic mode pulls toward start size

    BatchConfig cfg;
    cfg.num_trajectories   = (int)d->f_num_traj;
    cfg.initial_size       = (int)d->f_initial;
    cfg.target_generations = (long)d->f_target_gens;
    cfg.seed_base          = (unsigned int)time(NULL);
    cfg.base_params        = p;
    cfg.nbins              = (int)d->f_nbins;

    batch_init(&d->batch, cfg, 0);
    batch_start(&d->batch);
    d->has_batch = true;
    d->started = true;
}

// A labelled slider row; returns the (possibly updated) value via the pointer.
static void slider_row(float panel_x, float y, float w, const char *label,
                       const char *valtext, float *v, float lo, float hi) {
    DrawText(label, (int)panel_x + 12, (int)y + 2, 14, LIGHTGRAY);
    GuiSlider((Rectangle){panel_x + 130, y, w, 18}, NULL, valtext, v, lo, hi);
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

    // Ghost overlays: A. thaliana (real) only for now; uniform-model markers removed.
    RefMark ref_uniq[] = {{REF_REAL_UNIQUE_PER_KB, REAL, "real"}};
    RefMark ref_hors[] = {{REF_REAL_HORS_PER_KB, REAL, "real"}};
    RefMark ref_bs[]   = {{REF_REAL_BLOCK_SIZE_MED, REAL, "real med"}};
    RefMark ref_si[]   = {{REF_REAL_SIMILARITY_MED, REAL, "real med"}};

    struct { const char *title; HistSnap *s; RefMark *refs; int nrefs; } plots[6] = {
        {"UNIQUE REPEATS / kb (survivors)", &su, ref_uniq, 1},
        {"HORs / kb (survivors)",           &sh, ref_hors, 1},
        {"HOR BLOCK SIZE",                  &bs, ref_bs,   1},
        {"HOR BLOCK GAP",                   &bg, NULL,     0},
        {"HOR SIMILARITY",                  &si, ref_si,   1},
        {"COMPOSITE HOR METRIC",            &cp, NULL,     0},
    };
    Vector2 mouse = GetMousePosition();
    for (int i = 0; i < 6; i++) {
        draw_hist(cell[i], plots[i].title, plots[i].s, plots[i].refs, plots[i].nrefs,
                  d->plot_log_y[i], d->autoscale_x);
        // per-plot clickable Y-scale toggle (top-right of the cell)
        Rectangle tg = { cell[i].x + cell[i].width - 52, cell[i].y + 3, 46, 15 };
        bool hov = CheckCollisionPointRec(mouse, tg);
        DrawRectangleRec(tg, hov ? (Color){0,70,40,255} : (Color){25,35,25,255});
        DrawRectangleLinesEx(tg, 1, GRID);
        DrawText(d->plot_log_y[i] ? "Y:log" : "Y:lin", (int)tg.x + 5, (int)tg.y + 3, 9, BAR);
        if (hov && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) d->plot_log_y[i] = !d->plot_log_y[i];
    }

    // ----- status bar -----
    DrawRectangle(0, top - 2, plot_area_w + 4, 2, GRID);
    char status[256];
    if (d->has_batch) {
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
    GuiCheckBox((Rectangle){panel_x + 12, y, 20, 20}, "Unbounded", &d->unbounded); y += 34;

    bool busy = d->has_batch && (running || stopping);
    if (stopping) {
        GuiDisable();
        GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36}, "#133# Stopping...");
        GuiEnable();
    } else if (busy) {
        if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36}, "#133# Stop"))
            batch_request_stop(&d->batch);  // non-blocking; UI stays responsive
    } else {
        if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36}, "#131# Run batch"))
            launch(d);
    }
    y += 46;

    // Advanced (collapsible)
    if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 24},
                  d->show_advanced ? "#118# Advanced Options" : "#117# Advanced Options"))
        d->show_advanced = !d->show_advanced;
    y += 30;
    if (d->show_advanced) {
        slider_row(panel_x, y, sw, "Collapse <", TextFormat("%d", (int)d->f_collapse), &d->f_collapse, 0, 5000); y += 30;
        DrawText("Mutation", panel_x + 12, (int)y, 14, GRAY); y += 22;
        slider_row(panel_x, y, sw, "INDEL rate", TextFormat("%.2f", d->f_indel_rate), &d->f_indel_rate, 0.0f, 3.0f); y += 30;
        slider_row(panel_x, y, sw, "INDEL size", TextFormat("%.1f", d->f_indel_size), &d->f_indel_size, 1.0f, 100.0f); y += 30;
        slider_row(panel_x, y, sw, "SNP rate", TextFormat("%.2f", d->f_snp_rate), &d->f_snp_rate, 0.0f, 1.0f); y += 32;
        DrawText("Display", panel_x + 12, (int)y, 14, GRAY); y += 22;
        slider_row(panel_x, y, sw, "Bin resolution", TextFormat("%d", (int)d->f_nbins), &d->f_nbins, 20, 200); y += 26;
        DrawText("(applies on next Run)", panel_x + 14, (int)y, 10, (Color){120,100,60,255}); y += 18;
        GuiCheckBox((Rectangle){panel_x + 12, y, 20, 20}, "Autoscale X (fit data)", &d->autoscale_x); y += 26;
        DrawText("Click \"Y:log/lin\" on a plot to toggle", panel_x + 14, (int)y, 10, GRAY); y += 22;
    }
    y += 6;

    // legend
    DrawText("overlays:", panel_x + 12, (int)y, 11, GRAY); y += 16;
    DrawRectangle(panel_x + 12, (int)y + 2, 12, 8, REAL); DrawText("A. thaliana (report stand-in)", panel_x + 30, (int)y, 11, LIGHTGRAY); y += 16;
    DrawRectangle(panel_x + 12, (int)y + 2, 12, 8, MED);  DrawText("live median of data", panel_x + 30, (int)y, 11, LIGHTGRAY); y += 16;
    DrawText("see TODO.md: wire real data", panel_x + 12, (int)y, 10, (Color){120,100,60,255});
}

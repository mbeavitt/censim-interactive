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

#define DASH_MAXBINS 64

typedef struct {
    int   nbins, log_scale;
    float min, max;
    long  counts[DASH_MAXBINS];
    long  total, maxcount;
} HistSnap;

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
}

// Fraction [0,1] of a value across the histogram's x-range (log-aware).
static float frac_of(const HistSnap *s, float v) {
    float lo = s->min, hi = s->max, x = v;
    if (s->log_scale) {
        if (v <= 0) return -1.0f;
        lo = log10f(s->min); hi = log10f(s->max); x = log10f(v);
    }
    if (hi <= lo) return -1.0f;
    return (x - lo) / (hi - lo);
}

static const Color BG    = (Color){15, 20, 15, 255};
static const Color GRID  = (Color){0, 90, 40, 255};
static const Color BAR   = (Color){0, 200, 110, 255};
static const Color REAL  = (Color){255, 70, 70, 255};   // A. thaliana ghost
static const Color UNIF  = (Color){150, 150, 160, 255}; // uniform-model ghost

// Common chrome: frame, title, n=, x-scale tag. Returns the inner plot rect.
static Rectangle plot_frame(Rectangle b, const char *title, long n,
                            const char *xtag, const char *ytag) {
    DrawRectangleRec(b, BG);
    DrawRectangleLinesEx(b, 1, GRID);
    DrawText(title, (int)b.x + 6, (int)b.y + 4, 12, BAR);
    char nbuf[48];
    snprintf(nbuf, sizeof(nbuf), "n=%ld", n);
    int nw = MeasureText(nbuf, 10);
    DrawText(nbuf, (int)(b.x + b.width) - nw - 6, (int)b.y + 5, 10, GRID);
    if (xtag) DrawText(xtag, (int)(b.x + b.width) - MeasureText(xtag, 9) - 6, (int)(b.y + b.height) - 12, 9, GRID);
    if (ytag) DrawText(ytag, (int)b.x + 4, (int)b.y + 18, 9, GRID);
    return (Rectangle){ b.x + 4, b.y + 20, b.width - 8, b.height - 28 };
}

static void draw_hist(Rectangle b, const char *title, const HistSnap *s,
                      RefMark *refs, int nrefs, bool log_y) {
    const char *xtag = (s && s->log_scale) ? "log x" : "lin x";
    Rectangle p = plot_frame(b, title, s ? s->total : 0, xtag, log_y ? "log y" : "lin y");
    float px = p.x, py = p.y, pw = p.width, ph = p.height;

    if (!s || s->maxcount == 0) {
        DrawText("awaiting data", (int)(b.x + b.width/2 - 40), (int)(b.y + b.height/2), 10, GRID);
    } else {
        float denom = log_y ? log10f((float)s->maxcount + 1.0f) : (float)s->maxcount;
        float bw = pw / s->nbins;
        for (int i = 0; i < s->nbins; i++) {
            if (s->counts[i] == 0) continue;
            float num = log_y ? log10f((float)s->counts[i] + 1.0f) : (float)s->counts[i];
            float hh = ph * num / denom;
            DrawRectangle((int)(px + i*bw), (int)(py + ph - hh),
                          (int)(bw > 1 ? bw - 1 : 1), (int)hh, BAR);
        }
    }
    for (int r = 0; r < nrefs; r++) {
        float f = s ? frac_of(s, refs[r].value) : -1.0f;
        if (f < 0 || f > 1) continue;
        float x = px + f * pw;
        DrawLine((int)x, (int)py, (int)x, (int)(py + ph), refs[r].color);
        if (refs[r].label) DrawText(refs[r].label, (int)x + 2, (int)py + 1, 9, refs[r].color);
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
    d->log_y         = true;   // log-log on the count plots by default (matches paper)
    d->show_advanced = false;
    d->has_batch     = false;
    d->started       = false;
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
    int completed = 0, survived = 0, collapsed = 0, total = 0, workers = 0;
    bool running = false, complete = false;
    if (d->has_batch) {
        Batch *b = &d->batch;
        total = b->cfg.num_trajectories;
        workers = b->num_workers;
        pthread_mutex_lock(&b->lock);
        completed = b->completed; survived = b->survived; collapsed = b->collapsed_count;
        snap(&su, &b->h_unique_per_kb); snap(&sh, &b->h_hors_per_kb);
        snap(&bs, &b->h_block_size);    snap(&bg, &b->h_block_gap);
        snap(&si, &b->h_similarity);    snap(&dv, &b->h_diversity);
        snap(&cp, &b->h_composite);     snap(&cg, &b->h_collapse_gen);
        pthread_mutex_unlock(&b->lock);
        running = b->running && completed < total;
        complete = completed >= total;
        if (b->running && complete) batch_join(b);  // reap finished workers
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

    RefMark ref_uniq[] = {{REF_UNIFORM_UNIQUE_PER_KB, UNIF, "unif"}, {REF_REAL_UNIQUE_PER_KB, REAL, "real"}};
    RefMark ref_hors[] = {{REF_UNIFORM_HORS_PER_KB, UNIF, "unif"}, {REF_REAL_HORS_PER_KB, REAL, "real"}};
    RefMark ref_bs[]   = {{REF_REAL_BLOCK_SIZE_MED, REAL, "real med"}};
    RefMark ref_si[]   = {{REF_REAL_SIMILARITY_MED, REAL, "real med"}};

    bool ly = d->log_y;
    draw_hist(cell[0], "UNIQUE REPEATS / kb (survivors)", &su, ref_uniq, 2, false);  // paper: linear
    draw_hist(cell[1], "HORs / kb (survivors)",           &sh, ref_hors, 2, ly);     // paper: log x
    draw_hist(cell[2], "HOR BLOCK SIZE (log-log)",        &bs, ref_bs, 1, ly);       // paper: log-log
    draw_hist(cell[3], "HOR BLOCK GAP (log x)",           &bg, NULL, 0, ly);         // paper: log x
    draw_hist(cell[4], "HOR SIMILARITY",                  &si, ref_si, 1, false);    // paper: linear
    draw_hist(cell[5], "COMPOSITE HOR METRIC (log x)",    &cp, NULL, 0, ly);         // paper: log x

    // ----- status bar -----
    DrawRectangle(0, top - 2, plot_area_w + 4, 2, GRID);
    char status[256];
    if (d->has_batch) {
        float rate = total ? 100.0f * collapsed / total : 0.0f;
        snprintf(status, sizeof(status),
                 "%s  |  %d/%d done   survived %d   collapsed %d (%.0f%%)   |  %d workers",
                 running ? "RUNNING" : (complete ? "COMPLETE" : "READY"),
                 completed, total, survived, collapsed, rate, workers);
    } else {
        snprintf(status, sizeof(status), "Configure a batch and press Run >>");
    }
    DrawText(status, 8, 40, 14, running ? (Color){0,230,120,255} : LIGHTGRAY);

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
    GuiCheckBox((Rectangle){panel_x + 12, y, 20, 20}, "Unbounded (paper drift)", &d->unbounded); y += 34;

    bool busy = d->has_batch && running;
    if (busy) {
        if (GuiButton((Rectangle){panel_x + 12, y, panel_w - 24, 36}, "#133# Stop"))
            batch_stop(&d->batch);
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
        GuiCheckBox((Rectangle){panel_x + 12, y, 20, 20}, "Log Y (count plots)", &d->log_y); y += 30;
    }
    y += 6;

    // ghost legend
    DrawText("ghost overlays (report stand-ins):", panel_x + 12, (int)y, 11, GRAY); y += 16;
    DrawRectangle(panel_x + 12, (int)y + 2, 12, 8, REAL); DrawText("A. thaliana (real)", panel_x + 30, (int)y, 11, LIGHTGRAY); y += 16;
    DrawRectangle(panel_x + 12, (int)y + 2, 12, 8, UNIF); DrawText("uniform model", panel_x + 30, (int)y, 11, LIGHTGRAY); y += 16;
    DrawText("see TODO.md: wire real data", panel_x + 12, (int)y, 10, (Color){120,100,60,255});
}

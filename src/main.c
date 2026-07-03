#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#include <raylib.h>

// Forward declare for history functions
#include "simulation.h"
#include "stained_glass.h"   // SG_PAL_* ids used in the theme table below

// ============================================================================
// Theme + UI settings (colour schemes, UI scale, resizable panels)
// ============================================================================

typedef struct {
    const char *name;
    Color bg;         // app background
    Color panel;      // panel / box background
    Color header;     // panel header strip / stat box
    Color accent;     // accent lines + titles
    Color text;       // primary text
    Color text_dim;   // secondary text
    Color gui_base;   // raygui control fill
    Color gui_border; // raygui control border
    int   sg_pal;     // stained-glass palette id
} Theme;

static const Theme THEMES[] = {
    { "Midnight",  { 30, 30, 35,255},{ 25, 25, 30,255},{ 40, 40, 44,255},{  0,180,  0,255},{235,235,235,255},{150,150,160,255},{ 45, 45, 52,255},{ 90, 90,100,255}, SG_PAL_TURBO   },
    { "Amber"    , { 18, 14,  8,255},{ 26, 20, 10,255},{ 40, 32, 14,255},{255,176,  0,255},{255,214,140,255},{180,140, 70,255},{ 40, 32, 14,255},{120, 90, 30,255}, SG_PAL_MAGMA   },
    { "Ocean",     { 14, 22, 32,255},{ 18, 28, 40,255},{ 26, 40, 56,255},{  0,200,220,255},{210,230,240,255},{130,160,180,255},{ 26, 40, 56,255},{ 50, 90,120,255}, SG_PAL_ICE     },
    { "Grape",     { 24, 18, 32,255},{ 32, 24, 44,255},{ 44, 32, 60,255},{200,120,255,255},{230,215,245,255},{160,140,185,255},{ 40, 30, 56,255},{110, 80,150,255}, SG_PAL_VIRIDIS },
    { "Paper",     {235,232,224,255},{245,243,237,255},{220,216,206,255},{ 40, 90,160,255},{ 30, 30, 34,255},{ 90, 90,100,255},{225,222,214,255},{160,156,146,255}, SG_PAL_MONO    },
};
#define NUM_THEMES ((int)(sizeof(THEMES) / sizeof(THEMES[0])))

static Theme g_theme;              // active theme (copied from THEMES[g_theme_idx])
static int   g_theme_idx  = 0;
float        g_ui_scale   = 1.0f;  // non-static: shared with dashboard.c for text scaling
static int   g_panel_width;        // runtime, resizable + scaled (init in main)
static int   g_tile_size;          // runtime, scaled (init in main)

// Font-scaling wrappers so raw text tracks the same UI scale as raygui controls.
// All raw DrawText/MeasureText calls in this file are routed through these.
static int scaled_font(int size) {
    int s = (int)(size * g_ui_scale + 0.5f);
    return s < 1 ? 1 : s;
}
static void DrawTextS(const char *t, int x, int y, int size, Color c) {
    DrawText(t, x, y, scaled_font(size), c);
}
static int MeasureTextS(const char *t, int size) {
    return MeasureText(t, scaled_font(size));
}

// ============================================================================
// Stats History for Mission Control plots
// ============================================================================

#define HISTORY_MAX 1000  // Max data points to store

typedef struct {
    int generation[HISTORY_MAX];
    float diversity[HISTORY_MAX];
    int unique[HISTORY_MAX];
    int array_size[HISTORY_MAX];
    // Measured (realized) dup/del statistics over each sample interval -- these can
    // diverge from the theoretical values derived from the ratio because of the
    // size floor (0-size no-ops), placement skips and hard bounds.
    float m_dup_size[HISTORY_MAX];  // units per duplication event
    float m_del_size[HISTORY_MAX];  // units per deletion event
    float m_dup_rate[HISTORY_MAX];  // duplication events per generation
    float m_del_rate[HISTORY_MAX];  // deletion events per generation
    int count;
    int sample_interval;  // Sample every N generations
    int last_sampled_gen;
    // Baselines for per-interval deltas
    int  last_dup_count, last_del_count;
    long last_dup_units, last_del_units;
} StatsHistory;

static void history_init(StatsHistory *h) {
    h->count = 0;
    h->sample_interval = 100;
    h->last_sampled_gen = -1;
    h->last_dup_count = h->last_del_count = 0;
    h->last_dup_units = h->last_del_units = 0;
}

static void history_clear(StatsHistory *h) {
    h->count = 0;
    h->last_sampled_gen = -1;
    h->last_dup_count = h->last_del_count = 0;
    h->last_dup_units = h->last_del_units = 0;
}

static void history_record(StatsHistory *h, Simulation *sim) {
    // Only sample at intervals to avoid too many points
    if (sim->stats.generation - h->last_sampled_gen < h->sample_interval) return;

    int unique = sim_count_unique(sim);
    float diversity = (sim->array.num_units > 0)
        ? (float)unique / (float)sim->array.num_units : 0.0f;

    // Realized dup/del size & rate over the interval since the last sample.
    int base_gen = h->last_sampled_gen < 0 ? 0 : h->last_sampled_gen;
    int dgen = sim->stats.generation - base_gen; if (dgen < 1) dgen = 1;
    int  d_dup_ev = sim->stats.dup_count - h->last_dup_count;
    int  d_del_ev = sim->stats.del_count - h->last_del_count;
    long d_dup_u  = sim->stats.dup_units - h->last_dup_units;
    long d_del_u  = sim->stats.del_units - h->last_del_units;
    float m_dup_size = d_dup_ev > 0 ? (float)d_dup_u / d_dup_ev : 0.0f;
    float m_del_size = d_del_ev > 0 ? (float)d_del_u / d_del_ev : 0.0f;
    float m_dup_rate = (float)d_dup_ev / dgen;
    float m_del_rate = (float)d_del_ev / dgen;

    // Shift data if buffer full
    if (h->count >= HISTORY_MAX) {
        for (int i = 0; i < HISTORY_MAX - 1; i++) {
            h->generation[i] = h->generation[i + 1];
            h->diversity[i] = h->diversity[i + 1];
            h->unique[i] = h->unique[i + 1];
            h->array_size[i] = h->array_size[i + 1];
            h->m_dup_size[i] = h->m_dup_size[i + 1];
            h->m_del_size[i] = h->m_del_size[i + 1];
            h->m_dup_rate[i] = h->m_dup_rate[i + 1];
            h->m_del_rate[i] = h->m_del_rate[i + 1];
        }
        h->count = HISTORY_MAX - 1;
    }

    h->generation[h->count] = sim->stats.generation;
    h->diversity[h->count] = diversity;
    h->unique[h->count] = unique;
    h->array_size[h->count] = sim->array.num_units;
    h->m_dup_size[h->count] = m_dup_size;
    h->m_del_size[h->count] = m_del_size;
    h->m_dup_rate[h->count] = m_dup_rate;
    h->m_del_rate[h->count] = m_del_rate;
    h->count++;
    h->last_sampled_gen = sim->stats.generation;
    h->last_dup_count = sim->stats.dup_count;
    h->last_del_count = sim->stats.del_count;
    h->last_dup_units = sim->stats.dup_units;
    h->last_del_units = sim->stats.del_units;
}

// Draw a single retro-style plot
static void draw_plot(Rectangle bounds, const char *title, float *values, int count,
                      float min_val, float max_val, Color line_color, Color grid_color) {
    // Background
    DrawRectangleRec(bounds, (Color){15, 20, 15, 255});
    DrawRectangleLinesEx(bounds, 2, grid_color);

    // Title with glow effect
    DrawTextS(title, bounds.x + 8, bounds.y + 4, 14, grid_color);
    DrawTextS(title, bounds.x + 7, bounds.y + 3, 14, line_color);

    // Grid lines (retro scanline effect)
    for (int i = 1; i < 4; i++) {
        float y = bounds.y + (bounds.height * i / 4);
        DrawLine(bounds.x + 1, y, bounds.x + bounds.width - 1, y, (Color){grid_color.r, grid_color.g, grid_color.b, 60});
    }
    for (int i = 1; i < 6; i++) {
        float x = bounds.x + (bounds.width * i / 6);
        DrawLine(x, bounds.y + 20, x, bounds.y + bounds.height - 5, (Color){grid_color.r, grid_color.g, grid_color.b, 60});
    }

    if (count < 2) {
        DrawTextS("AWAITING DATA...", bounds.x + bounds.width/2 - 60, bounds.y + bounds.height/2 - 8, 12, grid_color);
        return;
    }

    // Plot area (with padding)
    float plot_x = bounds.x + 5;
    float plot_y = bounds.y + 22;
    float plot_w = bounds.width - 10;
    float plot_h = bounds.height - 30;

    // Auto-scale if min/max are equal
    if (max_val <= min_val) {
        min_val = values[0];
        max_val = values[0];
        for (int i = 1; i < count; i++) {
            if (values[i] < min_val) min_val = values[i];
            if (values[i] > max_val) max_val = values[i];
        }
        // Add padding
        float range = max_val - min_val;
        if (range < 0.001f) range = 1.0f;
        min_val -= range * 0.1f;
        max_val += range * 0.1f;
    }

    // Draw line with glow
    float x_step = plot_w / (count - 1);
    for (int i = 0; i < count - 1; i++) {
        float x1 = plot_x + i * x_step;
        float x2 = plot_x + (i + 1) * x_step;
        float y1 = plot_y + plot_h - ((values[i] - min_val) / (max_val - min_val)) * plot_h;
        float y2 = plot_y + plot_h - ((values[i + 1] - min_val) / (max_val - min_val)) * plot_h;

        // Glow (draw wider line behind)
        DrawLineEx((Vector2){x1, y1}, (Vector2){x2, y2}, 3.0f, (Color){line_color.r, line_color.g, line_color.b, 80});
        DrawLineEx((Vector2){x1, y1}, (Vector2){x2, y2}, 1.5f, line_color);
    }

    // Current value display
    char val_buf[32];
    if (max_val > 100) {
        snprintf(val_buf, sizeof(val_buf), "%.0f", values[count - 1]);
    } else {
        snprintf(val_buf, sizeof(val_buf), "%.3f", values[count - 1]);
    }
    int val_w = MeasureTextS(val_buf, 12);
    DrawTextS(val_buf, bounds.x + bounds.width - val_w - 8, bounds.y + 5, 12, line_color);
}

// Draw integer values (convert to float for plotting)
static void draw_plot_int(Rectangle bounds, const char *title, int *values, int count,
                          int min_val, int max_val, Color line_color, Color grid_color) {
    static float float_buf[HISTORY_MAX];
    for (int i = 0; i < count && i < HISTORY_MAX; i++) {
        float_buf[i] = (float)values[i];
    }
    draw_plot(bounds, title, float_buf, count, (float)min_val, (float)max_val, line_color, grid_color);
}

// Compact plot of a MEASURED series against its THEORETICAL value (dashed white
// reference line). Shared y-scale spans both so the divergence is legible.
static void draw_measure_plot(Rectangle b, const char *title, const float *vals, int count,
                              float theo, Color line) {
    DrawRectangleRec(b, (Color){15, 20, 15, 235});
    DrawRectangleLinesEx(b, 1, (Color){line.r, line.g, line.b, 110});
    DrawTextS(title, (int)b.x + 5, (int)b.y + 3, 10, line);

    float px = b.x + 4, py = b.y + 16, pw = b.width - 8, ph = b.height - 30;

    // y-range covers the measured samples and the theoretical line, padded.
    float mn = theo, mx = theo;
    for (int i = 0; i < count; i++) {
        if (vals[i] < mn) mn = vals[i];
        if (vals[i] > mx) mx = vals[i];
    }
    float rng = mx - mn;
    if (rng < 1e-6f) { mn -= 0.5f; mx += 0.5f; }
    else             { mn -= rng * 0.12f; mx += rng * 0.12f; }
    rng = mx - mn; if (rng <= 0) rng = 1;

    // theoretical reference (dashed)
    float ty = py + ph - ((theo - mn) / rng) * ph;
    for (int x = (int)px; x < (int)(px + pw); x += 6)
        DrawLine(x, (int)ty, x + 3, (int)ty, (Color){255, 255, 255, 90});

    // measured trace
    if (count >= 2) {
        float xs = pw / (count - 1);
        for (int i = 0; i < count - 1; i++) {
            float y1 = py + ph - ((vals[i]     - mn) / rng) * ph;
            float y2 = py + ph - ((vals[i + 1] - mn) / rng) * ph;
            DrawLineEx((Vector2){px + i * xs, y1}, (Vector2){px + (i + 1) * xs, y2}, 1.5f, line);
        }
    }

    // readout: current measured value (left) vs theoretical (right)
    float cur = count > 0 ? vals[count - 1] : 0.0f;
    DrawTextS(TextFormat("%.2f", cur), (int)b.x + 5, (int)(b.y + b.height) - 12, 9, line);
    const char *tt = TextFormat("th %.2f", theo);
    DrawTextS(tt, (int)(b.x + b.width) - MeasureTextS(tt, 9) - 4, (int)(b.y + b.height) - 12, 9,
             (Color){210, 210, 210, 200});
}

// 2x2 cluster of measured-vs-theoretical dup/del plots (single-trajectory view).
static void draw_measure_cluster(const StatsHistory *h, const Simulation *sim, int x, int y) {
    float r   = sim->params.dup_del_size_ratio; if (r <= 0) r = 1.0f;
    float sr  = sqrtf(r);
    float mu  = sim->params.indel_rate;
    float lam = sim->params.indel_size_lambda;
    // Theoretical (naive) values implied by the parameters.
    float th_dup_size = lam * sr, th_del_size = lam / sr;
    float th_dup_rate = mu / (1.0f + r), th_del_rate = mu * r / (1.0f + r);

    const Color DUP = (Color){0, 220, 255, 255};   // cyan
    const Color DEL = (Color){255, 180, 0, 255};   // amber
    int pw = 150, ph = 64, gap = 6;

    DrawTextS("MEASURED vs THEORY (dup/del)", x, y - 13, 10, (Color){150, 170, 150, 255});
    draw_measure_plot((Rectangle){x,              y,            pw, ph}, "DUP SIZE (u/event)",
                      h->m_dup_size, h->count, th_dup_size, DUP);
    draw_measure_plot((Rectangle){x + pw + gap,   y,            pw, ph}, "DEL SIZE (u/event)",
                      h->m_del_size, h->count, th_del_size, DEL);
    draw_measure_plot((Rectangle){x,              y + ph + gap, pw, ph}, "DUP RATE (ev/gen)",
                      h->m_dup_rate, h->count, th_dup_rate, DUP);
    draw_measure_plot((Rectangle){x + pw + gap,   y + ph + gap, pw, ph}, "DEL RATE (ev/gen)",
                      h->m_del_rate, h->count, th_del_rate, DEL);
}

// Draw the mission control panel
static void draw_mission_control(StatsHistory *h, int screen_width, int screen_height,
                                  int panel_height, bool minimized, int ctrl_panel_width) {
    int panel_y = screen_height - panel_height;
    int panel_width = screen_width - ctrl_panel_width;  // Stop at control panel

    // Panel background with border
    DrawRectangle(0, panel_y, panel_width, panel_height, g_theme.panel);
    DrawLine(0, panel_y, panel_width, panel_y, (Color){0, 180, 0, 200});
    DrawLine(0, panel_y + 1, panel_width, panel_y + 1, (Color){0, 100, 0, 150});

    if (minimized) return;

    // Title (centered in stats panel area)
    const char *mc_title = "[ STATISTICS ]";
    int title_w = MeasureTextS(mc_title, 14);
    DrawTextS(mc_title, panel_width / 2 - title_w / 2 + 1, panel_y + 6, 14, (Color){0, 60, 0, 255});
    DrawTextS(mc_title, panel_width / 2 - title_w / 2, panel_y + 5, 14, (Color){0, 200, 0, 255});

    // Three plots side by side
    int plot_margin = 15;
    int plot_top = panel_y + 28;
    int plot_height = panel_height - 38;
    int available_width = screen_width - ctrl_panel_width - plot_margin * 4;
    int plot_width = available_width / 3;

    // Diversity plot (green)
    Rectangle div_rect = {plot_margin, plot_top, plot_width, plot_height};
    draw_plot(div_rect, "DIVERSITY", h->diversity, h->count, 0.0f, 1.0f,
              (Color){0, 255, 100, 255}, (Color){0, 150, 60, 255});

    // Unique repeats plot (cyan)
    Rectangle uniq_rect = {plot_margin * 2 + plot_width, plot_top, plot_width, plot_height};
    draw_plot_int(uniq_rect, "UNIQUE SEQS", h->unique, h->count, 0, 0,
                  (Color){0, 220, 255, 255}, (Color){0, 120, 150, 255});

    // Array size plot (amber/orange)
    Rectangle size_rect = {plot_margin * 3 + plot_width * 2, plot_top, plot_width, plot_height};
    draw_plot_int(size_rect, "ARRAY SIZE", h->array_size, h->count, 0, 0,
                  (Color){255, 180, 0, 255}, (Color){150, 100, 0, 255});
}

// ============================================================================
// Resource path helpers (for app bundle support)
// ============================================================================

// Get the path to the visualize_umap executable/script
// In app bundle: uses CENSIM_RESOURCES env var pointing to bundled PyInstaller binary
// In dev mode: uses python3 with script relative to executable
// If output_path is NULL, uses --show mode to display interactively
static void get_umap_command(char *cmd, size_t cmd_size, const char *fasta_path,
                             const char *output_path, int grid_width) {
    const char *resources = getenv("CENSIM_RESOURCES");

    if (resources && strlen(resources) > 0) {
        // App bundle mode: use bundled PyInstaller executable
        if (output_path) {
            snprintf(cmd, cmd_size,
                "\"%s/visualize_umap/visualize_umap\" \"%s\" -o \"%s\" -w %d &",
                resources, fasta_path, output_path, grid_width);
        } else {
            snprintf(cmd, cmd_size,
                "\"%s/visualize_umap/visualize_umap\" \"%s\" --show -w %d &",
                resources, fasta_path, grid_width);
        }
    } else {
        // Development mode: use python3 with script
        if (output_path) {
            snprintf(cmd, cmd_size,
                "python3 \"%s/../../scripts/visualize_umap.py\" \"%s\" -o \"%s\" -w %d &",
                GetApplicationDirectory(), fasta_path, output_path, grid_width);
        } else {
            snprintf(cmd, cmd_size,
                "python3 \"%s/../../scripts/visualize_umap.py\" \"%s\" --show -w %d &",
                GetApplicationDirectory(), fasta_path, grid_width);
        }
    }
}

// Open a native "save file" dialog and place the chosen path in `out`.
// Returns true if a path was selected/derived, false if the user cancelled.
// macOS uses osascript; Linux prefers zenity (GNOME) then kdialog (KDE), and
// falls back to a default name in the working directory if neither is present.
static bool get_save_fasta_path(char *out, size_t out_size, int generation) {
    char default_name[64];
    snprintf(default_name, sizeof(default_name), "censim_gen%d.fasta", generation);

    char cmd[512];
#ifdef __APPLE__
    snprintf(cmd, sizeof(cmd),
        "osascript -e 'POSIX path of (choose file name with prompt \"Save FASTA as:\" default name \"%s\")' 2>/dev/null",
        default_name);
#else
    if (system("command -v zenity >/dev/null 2>&1") == 0) {
        snprintf(cmd, sizeof(cmd),
            "zenity --file-selection --save --confirm-overwrite "
            "--title=\"Save FASTA as\" --filename=\"%s\" 2>/dev/null",
            default_name);
    } else if (system("command -v kdialog >/dev/null 2>&1") == 0) {
        snprintf(cmd, sizeof(cmd),
            "kdialog --getsavefilename \".\" \"%s\" 2>/dev/null", default_name);
    } else {
        // No dialog utility available: write to the working directory.
        snprintf(out, out_size, "%s", default_name);
        return true;
    }
#endif

    FILE *pipe = popen(cmd, "r");
    if (!pipe) return false;

    char path[1024] = {0};
    bool ok = false;
    if (fgets(path, sizeof(path), pipe)) {
        path[strcspn(path, "\n")] = 0;
        if (strlen(path) > 0) {
            snprintf(out, out_size, "%s", path);
            ok = true;
        }
    }
    pclose(pipe);
    return ok;
}

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

#include "config.h"
#include "simulation.h"
#include "colorizer.h"
#include "dashboard.h"
#include "stained_glass.h"

// Apply a colour-scheme preset: copy its colours into g_theme and push the
// chrome colours into raygui + the stained-glass palette.
static void apply_theme(int idx, StainedGlass *sg) {
    if (idx < 0) idx = 0; if (idx >= NUM_THEMES) idx = NUM_THEMES - 1;
    g_theme_idx = idx;
    g_theme = THEMES[idx];

    GuiSetStyle(DEFAULT, BACKGROUND_COLOR,     ColorToInt(g_theme.panel));
    GuiSetStyle(DEFAULT, BASE_COLOR_NORMAL,    ColorToInt(g_theme.gui_base));
    GuiSetStyle(DEFAULT, BASE_COLOR_FOCUSED,   ColorToInt(g_theme.header));
    GuiSetStyle(DEFAULT, BASE_COLOR_PRESSED,   ColorToInt(g_theme.accent));
    GuiSetStyle(DEFAULT, BORDER_COLOR_NORMAL,  ColorToInt(g_theme.gui_border));
    GuiSetStyle(DEFAULT, BORDER_COLOR_FOCUSED, ColorToInt(g_theme.accent));
    GuiSetStyle(DEFAULT, TEXT_COLOR_NORMAL,    ColorToInt(g_theme.text));
    GuiSetStyle(DEFAULT, TEXT_COLOR_FOCUSED,   ColorToInt(g_theme.text));
    GuiSetStyle(DEFAULT, TEXT_COLOR_PRESSED,   ColorToInt(g_theme.bg));
    GuiSetStyle(DEFAULT, LINE_COLOR,           ColorToInt(g_theme.gui_border));

    // Sliders draw their value text to the RIGHT of the bar, on the panel
    // background -- not inside a filled control. The DEFAULT pressed/focused
    // text color (g_theme.bg) is meant for buttons whose base turns to accent,
    // so on a slider it makes the value vanish while you drag. Keep the slider
    // value legible in every state.
    GuiSetStyle(SLIDER, TEXT_COLOR_FOCUSED,    ColorToInt(g_theme.text));
    GuiSetStyle(SLIDER, TEXT_COLOR_PRESSED,    ColorToInt(g_theme.text));

    if (sg) sg_set_palette(sg, g_theme.sg_pal);
}

// Base (design) sizes to scale from, so scaling never compounds.
enum { UI_BASE_PANEL = PANEL_WIDTH, UI_BASE_TILE = TILE_SIZE };

// Apply a UI scale: raygui text size, array tile size, and (rescaled) panel width.
static void apply_ui_scale(float scale) {
    g_ui_scale = scale;
    GuiSetStyle(DEFAULT, TEXT_SIZE, (int)(16 * scale + 0.5f));
    g_tile_size   = (int)(UI_BASE_TILE * scale + 0.5f);   if (g_tile_size < 2) g_tile_size = 2;
    g_panel_width = (int)(UI_BASE_PANEL * scale + 0.5f);
}

// Thin vertical grip drawn on the control-panel splitter so it's discoverable.
static void draw_panel_grip(int x, int h, bool hot) {
    Color c = hot ? g_theme.accent : g_theme.gui_border;
    DrawRectangle(x - 1, 0, 2, h, c);
    int cy = h / 2;
    for (int i = -1; i <= 1; i++)
        DrawRectangle(x - 1, cy + i * 8, 2, 4, g_theme.accent);
}

// Options button + modal (UI scale, colour scheme). Uses function-local static
// state so it can be dropped into either view's draw path; drawn last = topmost.
static void draw_options_ui(int sw, int sh, StainedGlass *sg) {
    static const float SCALES[] = { 0.75f, 1.0f, 1.25f, 1.5f };
    static bool show = false;
    static int  scale_idx = 1;

    if (GuiButton((Rectangle){ (float)(sw - 116), (float)(sh - 38), 104, 26 }, "#142#Options")) show = true;
    if (!show) return;

    int mw = 440, mh = 250, margin = 20;
    Rectangle win = { (float)(sw - mw - margin), (float)(sh - mh - margin), (float)mw, (float)mh };
    DrawRectangle(0, 0, sw, sh, (Color){ 0, 0, 0, 130 });
    if (GuiWindowBox(win, "#142# Options")) { show = false; return; }

    int cx = (int)win.x + 16, cy = (int)win.y + 44;

    GuiLabel((Rectangle){ (float)cx, (float)cy, 200, 20 }, "UI scale");
    cy += 22;
    int sprev = scale_idx;
    GuiToggleGroup((Rectangle){ (float)cx, (float)cy, 64, 28 }, "0.75x;1.0x;1.25x;1.5x", &scale_idx);
    if (scale_idx != sprev) apply_ui_scale(SCALES[scale_idx]);
    cy += 50;

    GuiLabel((Rectangle){ (float)cx, (float)cy, 200, 20 }, "Colour scheme");
    cy += 22;
    char names[256] = "";
    for (int i = 0; i < NUM_THEMES; i++) {
        strcat(names, THEMES[i].name);
        if (i < NUM_THEMES - 1) strcat(names, ";");
    }
    int tprev = g_theme_idx;
    GuiToggleGroup((Rectangle){ (float)cx, (float)cy, 78, 28 }, names, &tprev);
    if (tprev != g_theme_idx) apply_theme(tprev, sg);
    cy += 46;

    DrawTextS("Tip: drag the panel's left edge or the self-identity",
             cx, cy, 11, g_theme.text_dim);
    DrawTextS("box corner to resize.", cx, cy + 14, 11, g_theme.text_dim);
}

// ============================================================================
// Grid rendering
// ============================================================================

static void draw_grid(Simulation *sim, Colorizer *colorizer, int offset_x, int offset_y,
                      int grid_width, int max_units) {
    int num_units = sim->array.num_units;
    if (num_units == 0) return;

    // Only draw tiles that fall inside the viewport. The grid is CPU/immediate-mode
    // (one hash + DrawRectangle per unit, every frame, even paused), so rendering the
    // rows scrolled off the bottom is pure waste -- for large arrays that is the
    // entire slowdown. Cap to the visible row count.
    if (max_units > 0 && num_units > max_units) num_units = max_units;

    for (int i = 0; i < num_units; i++) {
        int row = i / grid_width;
        int col = i % grid_width;

        Color c = colorizer_get_color(colorizer, sim->array.units[i]);

        int x = offset_x + col * g_tile_size;
        int y = offset_y + row * g_tile_size;

        DrawRectangle(x, y, g_tile_size - 1, g_tile_size - 1, c);
    }
}

// ============================================================================
// Stats panel
// ============================================================================

static void draw_stats(Simulation *sim, int unique, int x, int y, bool running) {
    // `unique` is passed in (cached, recomputed only when the array changes) so we
    // don't run the O(n) sim_count_unique on every frame.
    float diversity = (sim->array.num_units > 0)
        ? (float)unique / (float)sim->array.num_units
        : 0.0f;

    char buf[256];

    DrawRectangle(x, y, g_panel_width - 20, 215, g_theme.header);
    DrawRectangleLines(x, y, g_panel_width - 20, 215, LIGHTGRAY);

    int line = y + 15;
    int spacing = 22;

    DrawTextS("Statistics", x + 10, line, 20, WHITE);
    line += spacing + 10;

    snprintf(buf, sizeof(buf), "Generation: %d", sim->stats.generation);
    DrawTextS(buf, x + 10, line, 18, RAYWHITE);
    line += spacing;

    snprintf(buf, sizeof(buf), "Array size: %d", sim->array.num_units);
    DrawTextS(buf, x + 10, line, 18, RAYWHITE);
    line += spacing;

    snprintf(buf, sizeof(buf), "Unique seqs: %d", unique);
    DrawTextS(buf, x + 10, line, 18, RAYWHITE);
    line += spacing;

    snprintf(buf, sizeof(buf), "Diversity: %.4f", diversity);
    DrawTextS(buf, x + 10, line, 18, RAYWHITE);
    line += spacing + 10;

    DrawTextS("Mutations", x + 10, line, 16, GRAY);
    line += spacing - 4;

    snprintf(buf, sizeof(buf), "SNPs: %d  Dups: %d  Dels: %d",
             sim->stats.snp_count, sim->stats.dup_count, sim->stats.del_count);
    DrawTextS(buf, x + 10, line, 16, RAYWHITE);
    line += spacing;

    // Status indicator
    if (sim->stats.collapsed) {
        DrawTextS("COLLAPSED!", x + 10, line, 18, RED);
    } else if (running) {
        DrawTextS("RUNNING", x + 10, line, 18, GREEN);
    } else {
        DrawTextS("PAUSED", x + 10, line, 18, GRAY);
    }
}

// ============================================================================
// Main
// ============================================================================

int main(void) {
    // Initialize window at a safe small size first
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    InitWindow(800, 600, "Centromere Evolution Simulator");
    SetTargetFPS(60);

    // Initialise runtime UI-scale globals (g_panel_width / g_tile_size) before
    // anything reads them. Theme (which also sets the stained-glass palette) is
    // applied after the panel is created, below.
    apply_ui_scale(1.0f);

    // Get monitor dimensions
    int monitor = GetCurrentMonitor();
    int monitor_w = GetMonitorWidth(monitor);
    int monitor_h = GetMonitorHeight(monitor);
    printf("Monitor: %d, size: %dx%d\n", monitor, monitor_w, monitor_h);

    // Resize window to fit monitor and maximize
    SetWindowSize(monitor_w, monitor_h);
    MaximizeWindow();

    // Calculate grid based on monitor width
    int grid_width = (monitor_w - g_panel_width - 40) / g_tile_size;  // Extra margin
    if (grid_width < 10) grid_width = 10;
    printf("Grid width: %d tiles\n", grid_width);
    int last_screen_width = monitor_w;
    int last_panel_width = g_panel_width;
    int last_tile_size = g_tile_size;

    // Initialize simulation
    Simulation sim;
    sim_init(&sim, DEFAULT_INITIAL_SIZE, (unsigned int)time(NULL));
    // Snapshot the starting parameters so Reset can restore every slider to its
    // original value (sim_reset only rewinds the array/stats, not the params).
    SimParams initial_params = sim.params;

    // Initialize colorizer
    Colorizer colorizer;
    colorizer_init(&colorizer, 42);

    // Cached grid render. The coloured-square grid is expensive to build (a 712-d
    // colour projection per tile on cache misses) and was rebuilt every frame --
    // even while paused. Instead render it into a GPU texture once and re-render
    // only when the array actually changes (grid_dirty). Idle/paused frames then
    // cost a single textured-quad blit. cached_unique does the same for the O(n)
    // unique-sequence count that draw_stats showed every frame.
    RenderTexture2D grid_rt = LoadRenderTexture(monitor_w, monitor_h);
    SetTextureFilter(grid_rt.texture, TEXTURE_FILTER_POINT);
    bool grid_dirty = true;

    // Live self-identity ("stained glass") panel, bottom-left above the stats.
    StainedGlass stained_glass;
    sg_init(&stained_glass);
    int  cached_unique = 0;
    int  last_w_px = monitor_w, last_h_px = monitor_h;

    // UI state
    bool running = false;
    float gens_per_frame = 100.0f;

    // Resize-handle state
    bool  panel_dragging = false;    // dragging the control-panel splitter
    bool  sg_dragging = false;       // dragging the stained-glass box corner
    float sg_side = 0.0f;            // stained-glass box side in px (0 = init from height)
    int step_size = 10000;
    bool show_advanced = false;
    char step_size_text[16] = "10000";
    bool step_size_edit = false;
    double umap_start_time = 0.0;  // When UMAP was started (0 = not running)
    int refresh_counter = 0;
    float panel_scroll = 0.0f;  // Scroll offset for controls panel

    // Mission control state
    StatsHistory stats_history;
    history_init(&stats_history);
    bool mc_minimized = false;
    int mc_panel_height = 180;
    bool count_dist_edit = false;  // Dropdown state
    bool size_dist_edit = false;   // Dropdown state

    // App mode: 0 = single-trajectory view, 1 = multi-trajectory dashboard
    int app_mode = 0;
    Dashboard dash;
    dashboard_init(&dash);

    // raygui style + initial colour scheme (also sets the stained-glass palette)
    apply_theme(g_theme_idx, &stained_glass);

    printf("Centromere Evolution Simulator\n");
    printf("Controls:\n");
    printf("  Start/Stop: Toggle simulation\n");
    printf("  Step 1000: Advance 1000 generations\n");
    printf("  Reset: Restart simulation\n");
    printf("  Sliders: Adjust mutation rates\n");
    printf("\n");

    // Main loop
    while (!WindowShouldClose()) {
        // Toggle borderless fullscreen: Cmd+F on macOS, F11 on Windows/Linux.
        // ToggleBorderlessWindowed() resizes the window to the monitor without a
        // video-mode switch, so there's no resolution change or scaling artifacts
        // (raylib remembers the prior geometry to restore on toggle-off).
#ifdef __APPLE__
        bool toggle_fullscreen = (IsKeyDown(KEY_LEFT_SUPER) || IsKeyDown(KEY_RIGHT_SUPER)) && IsKeyPressed(KEY_F);
#else
        bool toggle_fullscreen = IsKeyPressed(KEY_F11);
#endif
        if (toggle_fullscreen) ToggleBorderlessWindowed();

        // Get current screen dimensions
        int screen_width = GetScreenWidth();
        int screen_height = GetScreenHeight();
        int panel_x = screen_width - g_panel_width;

        // Control-panel splitter: drag the panel's left edge to resize it.
        Vector2 mouse_pos = GetMousePosition();
        Rectangle splitter = { (float)panel_x - 4, 0, 8, (float)screen_height };
        bool over_splitter = CheckCollisionPointRec(mouse_pos, splitter);
        if (over_splitter && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) panel_dragging = true;
        if (panel_dragging) {
            g_panel_width = screen_width - (int)mouse_pos.x;
            if (g_panel_width < 300) g_panel_width = 300;
            if (g_panel_width > screen_width - 200) g_panel_width = screen_width - 200;
            panel_x = screen_width - g_panel_width;
            if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) panel_dragging = false;
        }
        if (over_splitter || panel_dragging) SetMouseCursor(MOUSE_CURSOR_RESIZE_EW);
        else SetMouseCursor(MOUSE_CURSOR_DEFAULT);

        // Recalculate grid width if the screen, panel width, or tile size changed.
        // (Panel width and tile size change via UI scale or the panel splitter.)
        if (abs(screen_width - last_screen_width) > 100 ||
            g_panel_width != last_panel_width || g_tile_size != last_tile_size) {
            grid_width = (screen_width - g_panel_width - 20) / g_tile_size;
            if (grid_width < 10) grid_width = 10;
            last_screen_width = screen_width;
            last_panel_width  = g_panel_width;
            last_tile_size    = g_tile_size;
            grid_dirty = true;
        }
        // Any size change alters the visible grid region -> rebuild the cache once.
        if (screen_width != last_w_px || screen_height != last_h_px) {
            grid_dirty = true; last_w_px = screen_width; last_h_px = screen_height;
        }

        // Update simulation (single-view only)
        if (app_mode == 0 && running && !sim.stats.collapsed) {
            sim_run(&sim, (int)gens_per_frame);
            grid_dirty = true;  // array changed -> grid + unique count are stale

            // Record stats for mission control
            history_record(&stats_history, &sim);

            // Refresh colorizer cache periodically
            refresh_counter++;
            if (refresh_counter >= 60) {  // Every second at 60 FPS
                refresh_counter = 0;
                // Don't clear cache - let it grow for performance
            }
        }

        // Re-render the grid into its cached GPU texture only when the array
        // changed. Skipped on idle/paused frames, so the per-tile colour projection
        // and the O(n) unique count stop running every frame -- the blit below is
        // all that remains on an idle frame.
        if (app_mode == 0 && grid_dirty) {
            cached_unique = sim_count_unique(&sim);
            int gmax = grid_width * ((screen_height - 10) / g_tile_size + 2);
            BeginTextureMode(grid_rt);
                ClearBackground(BLANK);
                draw_grid(&sim, &colorizer, 0, 0, grid_width, gmax);
            EndTextureMode();
            grid_dirty = false;
        }

        // Drawing
        BeginDrawing();
        ClearBackground(g_theme.bg);

        // Mode toggle button (top center, above the grid area)
        Rectangle toggle_rect = { (float)(screen_width - g_panel_width) / 2 - 80, 6, 160, 28 };

        // Dashboard mode: draw it and skip the single-trajectory view entirely
        if (app_mode == 1) {
            dashboard_update_draw(&dash, screen_width, screen_height, g_panel_width);
            if (GuiButton(toggle_rect, "#185# Single View")) app_mode = 0;
            DrawFPS(screen_width - 90, 6);
            draw_panel_grip(panel_x, screen_height, over_splitter || panel_dragging);
            draw_options_ui(screen_width, screen_height, &stained_glass);
            EndDrawing();
            continue;
        }

        // Blit the cached grid texture (rebuilt above only when dirty). The render
        // texture is stored bottom-up, so a negative source height flips it upright.
        DrawTextureRec(grid_rt.texture,
                       (Rectangle){ 0, 0, (float)grid_rt.texture.width, -(float)grid_rt.texture.height },
                       (Vector2){ 10, 10 }, WHITE);

        // Measured dup/del statistics vs theory (overlay, top-left of the grid)
        draw_measure_cluster(&stats_history, &sim, 12, 28);

        // Live self-identity panel, bottom-left, sitting just above the mission
        // control panel. Square, resizable via its top-right corner. ~30 fps.
        {
            int mc_h = mc_minimized ? 24 : mc_panel_height;
            float bottom = screen_height - mc_h - 12.0f;   // fixed bottom edge
            float max_side = bottom - 40.0f;               // leave a top margin
            if (sg_side <= 0.0f) sg_side = 0.45f * screen_height;   // first-frame default
            if (sg_side > max_side) sg_side = max_side;

            // Top-right corner drag handle.
            float bx = 12.0f, by = bottom - sg_side;
            Rectangle handle = { bx + sg_side - 16, by, 16, 16 };
            bool over_handle = CheckCollisionPointRec(mouse_pos, handle);
            if (over_handle && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) sg_dragging = true;
            if (sg_dragging) {
                sg_side = mouse_pos.x - bx;                 // keep square off the x drag
                if (sg_side < 150.0f) sg_side = 150.0f;
                if (sg_side > max_side) sg_side = max_side;
                by = bottom - sg_side;
                if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) sg_dragging = false;
                SetMouseCursor(MOUSE_CURSOR_RESIZE_NWSE);
            } else if (over_handle) {
                SetMouseCursor(MOUSE_CURSOR_RESIZE_NWSE);
            }

            Rectangle sg_box = { bx, by, sg_side, sg_side };
            sg_update_draw(&stained_glass, &sim, sg_box, 1.0 / 30.0,   // ~30 fps target
                           g_theme.accent, g_theme.panel, g_theme.text_dim);

            // Corner grip.
            Color grip = (over_handle || sg_dragging) ? g_theme.accent : g_theme.text_dim;
            for (int i = 0; i < 3; i++)
                DrawLine((int)(bx + sg_side) - 3 - i * 4, (int)by + 3,
                         (int)(bx + sg_side) - 3, (int)by + 3 + i * 4, grip);
        }

        // Draw panel background
        DrawRectangle(panel_x, 0, g_panel_width, screen_height, g_theme.panel);

        // Handle scrolling when mouse is over panel
        Rectangle panel_rect = {panel_x, 0, g_panel_width, screen_height};
        Vector2 mouse = GetMousePosition();

        // Calculate content height (approximate based on controls)
        int content_height = 830;  // Base height (incl. dup/del ratio row + 2-line readout)
        if (show_advanced) {
            content_height += 150;  // advanced block (freq derived from Dup/del ratio; no bias row/warning)
            if (sim.params.count_dist == DIST_NEGATIVE_BINOMIAL) content_height += 26;
            if (sim.params.size_dist == SIZE_POWER_LAW) content_height += 26;
        }

        // Max scroll is negative (scrolling down moves content up)
        float max_scroll = -(content_height - screen_height + 80);
        if (max_scroll > 0) max_scroll = 0;

        if (CheckCollisionPointRec(mouse, panel_rect)) {
            float wheel = GetMouseWheelMove();
            panel_scroll += wheel * 30.0f;
            // Clamp scroll
            if (panel_scroll > 0) panel_scroll = 0;
            if (panel_scroll < max_scroll) panel_scroll = max_scroll;
        }

        // Begin scissor mode for panel clipping
        BeginScissorMode(panel_x, 0, g_panel_width, screen_height);

        // Apply scroll offset to all controls
        int scroll_y = (int)panel_scroll;

        // Title (fixed, doesn't scroll)
        EndScissorMode();
        DrawRectangle(panel_x, 0, g_panel_width, 65, g_theme.header);
        DrawTextS("Controls", panel_x + 20, 20, 24, g_theme.text);
#ifdef __APPLE__
        DrawTextS("(Cmd+F fullscreen)", panel_x + 20, 48, 12, GRAY);
#else
        DrawTextS("(F11 fullscreen)", panel_x + 20, 48, 12, GRAY);
#endif
        BeginScissorMode(panel_x, 65, g_panel_width, screen_height - 65);

        // Buttons
        int btn_y = 75 + scroll_y;
        int btn_h = 40;
        int btn_spacing = 50;

        if (GuiButton((Rectangle){panel_x + 20, btn_y, 180, btn_h},
                      running ? "#132#Stop" : "#131#Start")) {
            running = !running;
        }
        if (GuiButton((Rectangle){panel_x + 210, btn_y, 180, btn_h}, "#72#Reset")) {
            running = false;
            sim_reset(&sim);
            sim.params = initial_params;   // restore all sliders to their originals
            gens_per_frame = 100.0f;
            step_size = 10000;
            snprintf(step_size_text, sizeof(step_size_text), "%d", step_size);
            colorizer_clear_cache(&colorizer);
            history_clear(&stats_history);
            grid_dirty = true;
        }
        btn_y += btn_spacing;

        if (GuiButton((Rectangle){panel_x + 20, btn_y, 180, btn_h}, TextFormat("#79#Step %d", step_size))) {
            sim_run(&sim, step_size);
            grid_dirty = true;
        }
        if (GuiButton((Rectangle){panel_x + 210, btn_y, 180, btn_h}, "#07#Export FASTA")) {
            char filepath[1024] = {0};
            if (get_save_fasta_path(filepath, sizeof(filepath), sim.stats.generation)) {
                FILE *f = fopen(filepath, "w");
                if (f) {
                    for (int i = 0; i < sim.array.num_units; i++) {
                        fprintf(f, ">repeat_%d\n%s\n", i + 1, sim.array.units[i]);
                    }
                    fclose(f);
                    printf("Exported %d repeats to %s\n", sim.array.num_units, filepath);
                } else {
                    fprintf(stderr, "Export FASTA: could not open '%s' for writing\n", filepath);
                }
            }
        }
        btn_y += btn_spacing + 20;

        // Sliders
        DrawTextS("Parameters", panel_x + 20, btn_y, 20, WHITE);
        btn_y += 30;

        int slider_w = 230;  // Leave room for value text on right
        int slider_h = 20;
        // Slider x-offset scales with the UI scale so the (scaled) labels to its
        // left don't overlap the sliders at larger scales.
        int label_w = (int)(100 * g_ui_scale + 0.5f);
        int row_h = 35;
        const char *hover_text = NULL;
        Rectangle hover_rect = {0};

        // INDEL rate
        Rectangle indel_row = {panel_x, btn_y - 5, g_panel_width, row_h};
        DrawTextS("INDEL rate:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%.2f", sim.params.indel_rate),
            &sim.params.indel_rate, 0.0f, 3.0f);
        if (CheckCollisionPointRec(mouse, indel_row)) {
            hover_text = "Expected INDELs (dup/del) per generation";
            hover_rect = indel_row;
        }
        btn_y += row_h;

        // Mean indel size (central / geometric-mean event size in repeat units;
        // the dup/del split is set by Dup/del ratio below).
        Rectangle size_row = {panel_x, btn_y - 5, g_panel_width, row_h};
        DrawTextS("Mean size:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%.1f", sim.params.indel_size_lambda),
            &sim.params.indel_size_lambda, 1.0f, 100.0f);
        if (CheckCollisionPointRec(mouse, size_row)) {
            hover_text = "Mean indel event size in repeat units (split into dup/del by Dup/del ratio)";
            hover_rect = size_row;
        }
        btn_y += row_h;

        // Dup/del ratio: a single dup:del SIZE ratio r (log slider centred at 1.0,
        // r in [0.001,1000]). Frequency is coupled to it -- dup_bias = 1/(1+r) -- so
        // that P(dup)*size_dup == P(del)*size_del and the mean array length is held
        // constant: bigger events are made proportionally rarer. The four derived
        // quantities (dup/del size and rate) are shown beneath.
        Rectangle ratio_row = {panel_x, btn_y - 5, g_panel_width, row_h};
        // Wrap "ratio" onto a second line so the label doesn't crowd the slider.
        DrawTextS("Dup/del", panel_x + 20, btn_y - 5, 16, LIGHTGRAY);
        DrawTextS("ratio:",  panel_x + 20, btn_y - 5 + scaled_font(16), 16, LIGHTGRAY);
        float size_ratio_e = log10f(sim.params.dup_del_size_ratio);
        const char *rfmt = sim.params.dup_del_size_ratio < 1.0f ? "%.3fx"
                         : sim.params.dup_del_size_ratio < 10.0f ? "%.2fx" : "%.0fx";
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat(rfmt, sim.params.dup_del_size_ratio),
            &size_ratio_e, -3.0f, 3.0f);
        sim.params.dup_del_size_ratio = powf(10.0f, size_ratio_e);
        sim.params.dup_bias = 1.0f / (1.0f + sim.params.dup_del_size_ratio);  // coupled freq
        if (CheckCollisionPointRec(mouse, ratio_row)) {
            hover_text = "Dup:del size ratio; frequency auto-balanced to hold array size (bigger = rarer)";
            hover_rect = ratio_row;
        }
        btn_y += row_h;

        // Four derived stats from (Mean size, INDEL rate, Dup/del ratio).
        {
            float r = sim.params.dup_del_size_ratio, sr = sqrtf(r);
            float pdup = 1.0f / (1.0f + r);
            DrawTextS(TextFormat("dup ~%.1f u  %.2f/gen", sim.params.indel_size_lambda * sr,
                                sim.params.indel_rate * pdup),
                     panel_x + label_w + 20, btn_y - 2, 12, (Color){130, 160, 130, 255});
            btn_y += 14;
            DrawTextS(TextFormat("del ~%.1f u  %.2f/gen", sim.params.indel_size_lambda / sr,
                                sim.params.indel_rate * (1.0f - pdup)),
                     panel_x + label_w + 20, btn_y - 2, 12, (Color){130, 160, 130, 255});
            btn_y += 18;
        }

        // SNP rate
        Rectangle snp_row = {panel_x, btn_y - 5, g_panel_width, row_h};
        DrawTextS("SNP rate:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%.2f", sim.params.snp_rate),
            &sim.params.snp_rate, 0.0f, 1.0f);
        if (CheckCollisionPointRec(mouse, snp_row)) {
            hover_text = "Expected point mutations per generation";
            hover_rect = snp_row;
        }
        btn_y += row_h;

        // Gens per frame
        Rectangle gpf_row = {panel_x, btn_y - 5, g_panel_width, row_h + 10};
        DrawTextS("Gens/frame:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%d", (int)gens_per_frame),
            &gens_per_frame, 10.0f, 100000.0f);
        if (CheckCollisionPointRec(mouse, gpf_row)) {
            hover_text = "Generations simulated per frame (speed control)";
            hover_rect = gpf_row;
        }
        btn_y += row_h + 10;

        // Target size
        Rectangle target_row = {panel_x, btn_y - 5, g_panel_width, row_h};
        DrawTextS("Target size:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        float target_f = (float)sim.params.target_size;
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%d", sim.params.target_size),
            &target_f, 1000.0f, 50000.0f);
        sim.params.target_size = (int)target_f;
        if (CheckCollisionPointRec(mouse, target_row)) {
            hover_text = "Target array size for elastic bounding";
            hover_rect = target_row;
        }
        btn_y += row_h;

        // Elasticity
        Rectangle elast_row = {panel_x, btn_y - 5, g_panel_width, row_h + 10};
        DrawTextS("Elasticity:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%.2f", sim.params.elasticity),
            &sim.params.elasticity, 0.0f, 1.0f);
        if (CheckCollisionPointRec(mouse, elast_row)) {
            hover_text = "Pull strength toward target size";
            hover_rect = elast_row;
        }
        btn_y += row_h + 10;

        // UMAP Visualization button
        // Check if UMAP signaled ready (file exists)
        bool umap_ready = (access("/tmp/censim_umap_ready", F_OK) == 0);
        if (umap_ready && umap_start_time > 0.0) {
            remove("/tmp/censim_umap_ready");
            umap_start_time = 0.0;  // Reset
        }
        bool umap_running = (umap_start_time > 0.0 && (GetTime() - umap_start_time) < 60.0);
        if (umap_running) {
            // Show running state
            GuiDisable();
            GuiButton((Rectangle){panel_x + 20, btn_y, 370, btn_h}, "#27#UMAP Running...");
            GuiEnable();
        } else if (GuiButton((Rectangle){panel_x + 20, btn_y, 370, btn_h}, "#27#UMAP Visualization (slow)")) {
            // Remove old signal file if exists
            remove("/tmp/censim_umap_ready");

            // Write temp FASTA
            char tempfasta[256];
            snprintf(tempfasta, sizeof(tempfasta), "/tmp/censim_temp_%d.fasta", sim.stats.generation);
            FILE *f = fopen(tempfasta, "w");
            if (f) {
                for (int i = 0; i < sim.array.num_units; i++) {
                    fprintf(f, ">repeat_%d\n%s\n", i + 1, sim.array.units[i]);
                }
                fclose(f);

                // Build and run visualization command (show mode)
                char script_cmd[2048];
                get_umap_command(script_cmd, sizeof(script_cmd), tempfasta, NULL, grid_width);
                printf("Running: %s\n", script_cmd);
                system(script_cmd);
                umap_start_time = GetTime();
            }
        }
        btn_y += btn_spacing;

        // Stats panel
        draw_stats(&sim, cached_unique, panel_x + 10, btn_y, running);
        btn_y += 225;

        // Advanced options (collapsible)
        if (GuiButton((Rectangle){panel_x + 20, btn_y, 370, 25},
                      show_advanced ? "#120#Advanced Options" : "#119#Advanced Options")) {
            show_advanced = !show_advanced;
        }
        btn_y += 30;

        if (show_advanced) {
            // Calculate dynamic height based on visible parameter sliders
            int adv_height = 145;
            if (sim.params.count_dist == DIST_NEGATIVE_BINOMIAL) adv_height += 26;
            if (sim.params.size_dist == SIZE_POWER_LAW) adv_height += 26;

            DrawRectangle(panel_x + 10, btn_y, g_panel_width - 20, adv_height, (Color){40, 40, 40, 200});

            int adv_y = btn_y + 10;

            // Step size
            Rectangle step_row = {panel_x + 20, adv_y - 2, 220, 24};
            DrawTextS("Step size:", panel_x + 20, adv_y, 16, LIGHTGRAY);
            if (GuiTextBox((Rectangle){panel_x + 120, adv_y - 3, 100, 24},
                          step_size_text, 16, step_size_edit)) {
                step_size_edit = !step_size_edit;
            }
            if (!step_size_edit) {
                int val = atoi(step_size_text);
                if (val > 0) step_size = val;
            }
            if (CheckCollisionPointRec(mouse, step_row)) {
                hover_text = "Generations per 'Step N' button click";
            }
            adv_y += 40;

            // (Dup/del frequency is derived from the Dup/del ratio above to keep the
            // mean array length constant, so there's no separate bias control here.)

            // Hard bounds checkbox
            Rectangle bounds_row = {panel_x + 20, adv_y, 300, 20};
            GuiCheckBox((Rectangle){panel_x + 20, adv_y, 20, 20}, "Hard bounds (min/max)", &sim.params.bounding_enabled);
            if (CheckCollisionPointRec(mouse, bounds_row)) {
                hover_text = "Enforce min/max array size limits";
            }
            adv_y += 35;

            // Events dropdown row
            int events_y = adv_y;
            Rectangle events_row = {panel_x + 20, events_y, 250, 24};
            DrawTextS("Events:", panel_x + 20, events_y + 2, 16, LIGHTGRAY);
            if (CheckCollisionPointRec(mouse, events_row) && !count_dist_edit && !size_dist_edit) {
                hover_text = "Distribution for mutation event counts per generation";
            }
            adv_y += 28;

            // Dispersion slider row (only if NB selected)
            int dispersion_y = adv_y;
            if (sim.params.count_dist == DIST_NEGATIVE_BINOMIAL) {
                Rectangle disp_row = {panel_x + 40, dispersion_y, 300, 20};
                DrawTextS("dispersion:", panel_x + 40, dispersion_y + 2, 14, GRAY);
                GuiSlider((Rectangle){panel_x + 140, dispersion_y, 180, 18}, NULL, NULL,
                          &sim.params.nb_dispersion, 0.1f, 5.0f);
                DrawTextS(TextFormat("%.1f", sim.params.nb_dispersion), panel_x + 330, dispersion_y + 2, 12, WHITE);
                if (CheckCollisionPointRec(mouse, disp_row)) {
                    hover_text = "Overdispersion (k): lower = more variance, higher = more Poisson-like";
                }
                adv_y += 26;
            }

            // Sizes dropdown row
            int sizes_y = adv_y;
            Rectangle sizes_row = {panel_x + 20, sizes_y, 250, 24};
            DrawTextS("Sizes:", panel_x + 20, sizes_y + 2, 16, LIGHTGRAY);
            if (CheckCollisionPointRec(mouse, sizes_row) && !count_dist_edit && !size_dist_edit) {
                hover_text = "Distribution for mutation sizes (dup/del length)";
            }
            adv_y += 28;

            // Alpha slider row (only if power law selected)
            if (sim.params.size_dist == SIZE_POWER_LAW) {
                Rectangle alpha_row = {panel_x + 40, adv_y, 300, 20};
                DrawTextS("alpha:", panel_x + 40, adv_y + 2, 14, GRAY);
                GuiSlider((Rectangle){panel_x + 140, adv_y, 180, 18}, NULL, NULL,
                          &sim.params.power_law_alpha, 1.5f, 4.0f);
                DrawTextS(TextFormat("%.1f", sim.params.power_law_alpha), panel_x + 330, adv_y + 2, 12, WHITE);
                if (CheckCollisionPointRec(mouse, alpha_row)) {
                    hover_text = "Tail heaviness: lower = heavier tail (more large events)";
                }
            }

            // Draw dropdowns last - lock the other when one is open to prevent click-through
            int count_dist = (int)sim.params.count_dist;
            int size_dist = (int)sim.params.size_dist;

            if (count_dist_edit) {
                // Events is open - draw sizes as disabled, then events on top
                GuiDisable();
                GuiDropdownBox((Rectangle){panel_x + 100, sizes_y, 160, 24},
                               "Poisson;Geometric;Power Law", &size_dist, false);
                GuiEnable();
                if (GuiDropdownBox((Rectangle){panel_x + 100, events_y, 160, 24},
                                   "Poisson;Negative Binomial", &count_dist, count_dist_edit)) {
                    count_dist_edit = !count_dist_edit;
                }
            } else if (size_dist_edit) {
                // Sizes is open - draw events as disabled, then sizes on top
                GuiDisable();
                GuiDropdownBox((Rectangle){panel_x + 100, events_y, 160, 24},
                               "Poisson;Negative Binomial", &count_dist, false);
                GuiEnable();
                if (GuiDropdownBox((Rectangle){panel_x + 100, sizes_y, 160, 24},
                                   "Poisson;Geometric;Power Law", &size_dist, size_dist_edit)) {
                    size_dist_edit = !size_dist_edit;
                }
            } else {
                // Neither open - draw both normally
                if (GuiDropdownBox((Rectangle){panel_x + 100, events_y, 160, 24},
                                   "Poisson;Negative Binomial", &count_dist, count_dist_edit)) {
                    count_dist_edit = !count_dist_edit;
                }
                if (GuiDropdownBox((Rectangle){panel_x + 100, sizes_y, 160, 24},
                                   "Poisson;Geometric;Power Law", &size_dist, size_dist_edit)) {
                    size_dist_edit = !size_dist_edit;
                }
            }

            sim.params.count_dist = (CountDistribution)count_dist;
            sim.params.size_dist = (SizeDistribution)size_dist;

            btn_y += adv_height + 5;
        }

        // End scissor mode before drawing overlays
        EndScissorMode();

        // FPS counter
        DrawFPS(screen_width - 100, 10);

        // Draw hover tooltip
        if (hover_text) {
            int text_width = MeasureTextS(hover_text, 14);
            int tip_x = (int)mouse.x + 15;
            int tip_y = (int)mouse.y - 25;
            if (tip_x + text_width + 10 > screen_width) tip_x = screen_width - text_width - 15;
            if (tip_y < 5) tip_y = (int)mouse.y + 20;
            DrawRectangle(tip_x - 5, tip_y - 3, text_width + 10, 20, (Color){50, 50, 55, 240});
            DrawRectangleLines(tip_x - 5, tip_y - 3, text_width + 10, 20, GRAY);
            DrawTextS(hover_text, tip_x, tip_y, 14, WHITE);
        }

        // Statistics panel
        int mc_height = mc_minimized ? 24 : mc_panel_height;
        draw_mission_control(&stats_history, screen_width, screen_height, mc_height, mc_minimized, g_panel_width);

        // Minimize/maximize toggle button (far left)
        int mc_btn_x = 10;
        int mc_btn_y = screen_height - mc_height + 2;
        const char *mc_btn_text = mc_minimized ? "^ STATS" : "v STATS";
        if (GuiButton((Rectangle){mc_btn_x, mc_btn_y, 90, 20}, mc_btn_text)) {
            mc_minimized = !mc_minimized;
        }

        // Mode toggle to enter the multi-trajectory dashboard (drawn on top)
        if (GuiButton(toggle_rect, "#191# Dashboard")) app_mode = 1;

        draw_panel_grip(panel_x, screen_height, over_splitter || panel_dragging);
        draw_options_ui(screen_width, screen_height, &stained_glass);
        EndDrawing();
    }

    // Cleanup
    dashboard_free(&dash);
    sim_free(&sim);
    colorizer_free(&colorizer);
    sg_free(&stained_glass);
    UnloadRenderTexture(grid_rt);
    CloseWindow();

    return 0;
}

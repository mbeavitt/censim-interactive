#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#include <raylib.h>

// Forward declare for history functions
#include "simulation.h"

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
    DrawText(title, bounds.x + 8, bounds.y + 4, 14, grid_color);
    DrawText(title, bounds.x + 7, bounds.y + 3, 14, line_color);

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
        DrawText("AWAITING DATA...", bounds.x + bounds.width/2 - 60, bounds.y + bounds.height/2 - 8, 12, grid_color);
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
    int val_w = MeasureText(val_buf, 12);
    DrawText(val_buf, bounds.x + bounds.width - val_w - 8, bounds.y + 5, 12, line_color);
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
    DrawText(title, (int)b.x + 5, (int)b.y + 3, 10, line);

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
    DrawText(TextFormat("%.2f", cur), (int)b.x + 5, (int)(b.y + b.height) - 12, 9, line);
    const char *tt = TextFormat("th %.2f", theo);
    DrawText(tt, (int)(b.x + b.width) - MeasureText(tt, 9) - 4, (int)(b.y + b.height) - 12, 9,
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

    DrawText("MEASURED vs THEORY (dup/del)", x, y - 13, 10, (Color){150, 170, 150, 255});
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
    DrawRectangle(0, panel_y, panel_width, panel_height, (Color){20, 25, 20, 245});
    DrawLine(0, panel_y, panel_width, panel_y, (Color){0, 180, 0, 200});
    DrawLine(0, panel_y + 1, panel_width, panel_y + 1, (Color){0, 100, 0, 150});

    if (minimized) return;

    // Title (centered in stats panel area)
    const char *mc_title = "[ STATISTICS ]";
    int title_w = MeasureText(mc_title, 14);
    DrawText(mc_title, panel_width / 2 - title_w / 2 + 1, panel_y + 6, 14, (Color){0, 60, 0, 255});
    DrawText(mc_title, panel_width / 2 - title_w / 2, panel_y + 5, 14, (Color){0, 200, 0, 255});

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

        int x = offset_x + col * TILE_SIZE;
        int y = offset_y + row * TILE_SIZE;

        DrawRectangle(x, y, TILE_SIZE - 1, TILE_SIZE - 1, c);
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

    DrawRectangle(x, y, PANEL_WIDTH - 20, 215, (Color){40, 40, 40, 220});
    DrawRectangleLines(x, y, PANEL_WIDTH - 20, 215, LIGHTGRAY);

    int line = y + 15;
    int spacing = 22;

    DrawText("Statistics", x + 10, line, 20, WHITE);
    line += spacing + 10;

    snprintf(buf, sizeof(buf), "Generation: %d", sim->stats.generation);
    DrawText(buf, x + 10, line, 18, RAYWHITE);
    line += spacing;

    snprintf(buf, sizeof(buf), "Array size: %d", sim->array.num_units);
    DrawText(buf, x + 10, line, 18, RAYWHITE);
    line += spacing;

    snprintf(buf, sizeof(buf), "Unique seqs: %d", unique);
    DrawText(buf, x + 10, line, 18, RAYWHITE);
    line += spacing;

    snprintf(buf, sizeof(buf), "Diversity: %.4f", diversity);
    DrawText(buf, x + 10, line, 18, RAYWHITE);
    line += spacing + 10;

    DrawText("Mutations", x + 10, line, 16, GRAY);
    line += spacing - 4;

    snprintf(buf, sizeof(buf), "SNPs: %d  Dups: %d  Dels: %d",
             sim->stats.snp_count, sim->stats.dup_count, sim->stats.del_count);
    DrawText(buf, x + 10, line, 16, RAYWHITE);
    line += spacing;

    // Status indicator
    if (sim->stats.collapsed) {
        DrawText("COLLAPSED!", x + 10, line, 18, RED);
    } else if (running) {
        DrawText("RUNNING", x + 10, line, 18, GREEN);
    } else {
        DrawText("PAUSED", x + 10, line, 18, GRAY);
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

    // Get monitor dimensions
    int monitor = GetCurrentMonitor();
    int monitor_w = GetMonitorWidth(monitor);
    int monitor_h = GetMonitorHeight(monitor);
    printf("Monitor: %d, size: %dx%d\n", monitor, monitor_w, monitor_h);

    // Resize window to fit monitor and maximize
    SetWindowSize(monitor_w, monitor_h);
    MaximizeWindow();

    // Calculate grid based on monitor width
    int grid_width = (monitor_w - PANEL_WIDTH - 40) / TILE_SIZE;  // Extra margin
    if (grid_width < 10) grid_width = 10;
    printf("Grid width: %d tiles\n", grid_width);
    int last_screen_width = monitor_w;

    // Initialize simulation
    Simulation sim;
    sim_init(&sim, DEFAULT_INITIAL_SIZE, (unsigned int)time(NULL));

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
    int  cached_unique = 0;
    int  last_w_px = monitor_w, last_h_px = monitor_h;

    // UI state
    bool running = false;
    float gens_per_frame = 100.0f;
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

    // raygui style
    GuiSetStyle(DEFAULT, TEXT_SIZE, 16);

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
        int panel_x = screen_width - PANEL_WIDTH;

        // Recalculate grid width if screen size changed significantly (>100px)
        if (abs(screen_width - last_screen_width) > 100) {
            grid_width = (screen_width - PANEL_WIDTH - 20) / TILE_SIZE;
            if (grid_width < 10) grid_width = 10;
            last_screen_width = screen_width;
            printf("Grid width changed: %d tiles\n", grid_width);
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
            int gmax = grid_width * ((screen_height - 10) / TILE_SIZE + 2);
            BeginTextureMode(grid_rt);
                ClearBackground(BLANK);
                draw_grid(&sim, &colorizer, 0, 0, grid_width, gmax);
            EndTextureMode();
            grid_dirty = false;
        }

        // Drawing
        BeginDrawing();
        ClearBackground((Color){30, 30, 35, 255});

        // Mode toggle button (top center, above the grid area)
        Rectangle toggle_rect = { (float)(screen_width - PANEL_WIDTH) / 2 - 80, 6, 160, 28 };

        // Dashboard mode: draw it and skip the single-trajectory view entirely
        if (app_mode == 1) {
            dashboard_update_draw(&dash, screen_width, screen_height, PANEL_WIDTH);
            if (GuiButton(toggle_rect, "#185# Single View")) app_mode = 0;
            DrawFPS(screen_width - 90, 6);
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

        // Draw panel background
        DrawRectangle(panel_x, 0, PANEL_WIDTH, screen_height, (Color){25, 25, 30, 255});

        // Handle scrolling when mouse is over panel
        Rectangle panel_rect = {panel_x, 0, PANEL_WIDTH, screen_height};
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
        BeginScissorMode(panel_x, 0, PANEL_WIDTH, screen_height);

        // Apply scroll offset to all controls
        int scroll_y = (int)panel_scroll;

        // Title (fixed, doesn't scroll)
        EndScissorMode();
        DrawRectangle(panel_x, 0, PANEL_WIDTH, 65, (Color){25, 25, 30, 255});
        DrawText("Controls", panel_x + 20, 20, 24, WHITE);
#ifdef __APPLE__
        DrawText("(Cmd+F fullscreen)", panel_x + 20, 48, 12, GRAY);
#else
        DrawText("(F11 fullscreen)", panel_x + 20, 48, 12, GRAY);
#endif
        BeginScissorMode(panel_x, 65, PANEL_WIDTH, screen_height - 65);

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
        DrawText("Parameters", panel_x + 20, btn_y, 20, WHITE);
        btn_y += 30;

        int slider_w = 230;  // Leave room for value text on right
        int slider_h = 20;
        int label_w = 100;
        int row_h = 35;
        const char *hover_text = NULL;
        Rectangle hover_rect = {0};

        // INDEL rate
        Rectangle indel_row = {panel_x, btn_y - 5, PANEL_WIDTH, row_h};
        DrawText("INDEL rate:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
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
        Rectangle size_row = {panel_x, btn_y - 5, PANEL_WIDTH, row_h};
        DrawText("Mean size:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
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
        Rectangle ratio_row = {panel_x, btn_y - 5, PANEL_WIDTH, row_h};
        DrawText("Dup/del ratio:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
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
            DrawText(TextFormat("dup ~%.1f u  %.2f/gen", sim.params.indel_size_lambda * sr,
                                sim.params.indel_rate * pdup),
                     panel_x + label_w + 20, btn_y - 2, 12, (Color){130, 160, 130, 255});
            btn_y += 14;
            DrawText(TextFormat("del ~%.1f u  %.2f/gen", sim.params.indel_size_lambda / sr,
                                sim.params.indel_rate * (1.0f - pdup)),
                     panel_x + label_w + 20, btn_y - 2, 12, (Color){130, 160, 130, 255});
            btn_y += 18;
        }

        // SNP rate
        Rectangle snp_row = {panel_x, btn_y - 5, PANEL_WIDTH, row_h};
        DrawText("SNP rate:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
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
        Rectangle gpf_row = {panel_x, btn_y - 5, PANEL_WIDTH, row_h + 10};
        DrawText("Gens/frame:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
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
        Rectangle target_row = {panel_x, btn_y - 5, PANEL_WIDTH, row_h};
        DrawText("Target size:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
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
        Rectangle elast_row = {panel_x, btn_y - 5, PANEL_WIDTH, row_h + 10};
        DrawText("Elasticity:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
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

            DrawRectangle(panel_x + 10, btn_y, PANEL_WIDTH - 20, adv_height, (Color){40, 40, 40, 200});

            int adv_y = btn_y + 10;

            // Step size
            Rectangle step_row = {panel_x + 20, adv_y - 2, 220, 24};
            DrawText("Step size:", panel_x + 20, adv_y, 16, LIGHTGRAY);
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
            DrawText("Events:", panel_x + 20, events_y + 2, 16, LIGHTGRAY);
            if (CheckCollisionPointRec(mouse, events_row) && !count_dist_edit && !size_dist_edit) {
                hover_text = "Distribution for mutation event counts per generation";
            }
            adv_y += 28;

            // Dispersion slider row (only if NB selected)
            int dispersion_y = adv_y;
            if (sim.params.count_dist == DIST_NEGATIVE_BINOMIAL) {
                Rectangle disp_row = {panel_x + 40, dispersion_y, 300, 20};
                DrawText("dispersion:", panel_x + 40, dispersion_y + 2, 14, GRAY);
                GuiSlider((Rectangle){panel_x + 140, dispersion_y, 180, 18}, NULL, NULL,
                          &sim.params.nb_dispersion, 0.1f, 5.0f);
                DrawText(TextFormat("%.1f", sim.params.nb_dispersion), panel_x + 330, dispersion_y + 2, 12, WHITE);
                if (CheckCollisionPointRec(mouse, disp_row)) {
                    hover_text = "Overdispersion (k): lower = more variance, higher = more Poisson-like";
                }
                adv_y += 26;
            }

            // Sizes dropdown row
            int sizes_y = adv_y;
            Rectangle sizes_row = {panel_x + 20, sizes_y, 250, 24};
            DrawText("Sizes:", panel_x + 20, sizes_y + 2, 16, LIGHTGRAY);
            if (CheckCollisionPointRec(mouse, sizes_row) && !count_dist_edit && !size_dist_edit) {
                hover_text = "Distribution for mutation sizes (dup/del length)";
            }
            adv_y += 28;

            // Alpha slider row (only if power law selected)
            if (sim.params.size_dist == SIZE_POWER_LAW) {
                Rectangle alpha_row = {panel_x + 40, adv_y, 300, 20};
                DrawText("alpha:", panel_x + 40, adv_y + 2, 14, GRAY);
                GuiSlider((Rectangle){panel_x + 140, adv_y, 180, 18}, NULL, NULL,
                          &sim.params.power_law_alpha, 1.5f, 4.0f);
                DrawText(TextFormat("%.1f", sim.params.power_law_alpha), panel_x + 330, adv_y + 2, 12, WHITE);
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
            int text_width = MeasureText(hover_text, 14);
            int tip_x = (int)mouse.x + 15;
            int tip_y = (int)mouse.y - 25;
            if (tip_x + text_width + 10 > screen_width) tip_x = screen_width - text_width - 15;
            if (tip_y < 5) tip_y = (int)mouse.y + 20;
            DrawRectangle(tip_x - 5, tip_y - 3, text_width + 10, 20, (Color){50, 50, 55, 240});
            DrawRectangleLines(tip_x - 5, tip_y - 3, text_width + 10, 20, GRAY);
            DrawText(hover_text, tip_x, tip_y, 14, WHITE);
        }

        // Statistics panel
        int mc_height = mc_minimized ? 24 : mc_panel_height;
        draw_mission_control(&stats_history, screen_width, screen_height, mc_height, mc_minimized, PANEL_WIDTH);

        // Minimize/maximize toggle button (far left)
        int mc_btn_x = 10;
        int mc_btn_y = screen_height - mc_height + 2;
        const char *mc_btn_text = mc_minimized ? "^ STATS" : "v STATS";
        if (GuiButton((Rectangle){mc_btn_x, mc_btn_y, 90, 20}, mc_btn_text)) {
            mc_minimized = !mc_minimized;
        }

        // Mode toggle to enter the multi-trajectory dashboard (drawn on top)
        if (GuiButton(toggle_rect, "#191# Dashboard")) app_mode = 1;

        EndDrawing();
    }

    // Cleanup
    dashboard_free(&dash);
    sim_free(&sim);
    colorizer_free(&colorizer);
    UnloadRenderTexture(grid_rt);
    CloseWindow();

    return 0;
}

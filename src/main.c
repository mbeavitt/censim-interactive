#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <raylib.h>

// ============================================================================
// Resource path helpers (for app bundle support)
// ============================================================================

// Get the path to the visualize_umap executable/script
// In app bundle: uses CENSIM_RESOURCES env var pointing to bundled PyInstaller binary
// In dev mode: uses python3 with script relative to executable
static void get_umap_command(char *cmd, size_t cmd_size, const char *fasta_path,
                             const char *output_path, int grid_width) {
    const char *resources = getenv("CENSIM_RESOURCES");

    if (resources && strlen(resources) > 0) {
        // App bundle mode: use bundled PyInstaller executable
        snprintf(cmd, cmd_size,
            "\"%s/visualize_umap/visualize_umap\" \"%s\" -o \"%s\" -w %d &",
            resources, fasta_path, output_path, grid_width);
    } else {
        // Development mode: use python3 with script
        snprintf(cmd, cmd_size,
            "python3 \"%s/../../scripts/visualize_umap.py\" \"%s\" -o \"%s\" -w %d &",
            GetApplicationDirectory(), fasta_path, output_path, grid_width);
    }
}

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

#include "config.h"
#include "simulation.h"
#include "colorizer.h"

// ============================================================================
// Grid rendering
// ============================================================================

static void draw_grid(Simulation *sim, Colorizer *colorizer, int offset_x, int offset_y, int grid_width) {
    int num_units = sim->array.num_units;
    if (num_units == 0) return;

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

static void draw_stats(Simulation *sim, int x, int y, bool running) {
    int unique = sim_count_unique(sim);
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
    sim_init(&sim, DEFAULT_INITIAL_SIZE);

    // Initialize colorizer
    Colorizer colorizer;
    colorizer_init(&colorizer, 42);

    // UI state
    bool running = false;
    float gens_per_frame = 100.0f;
    int step_size = 10000;
    bool show_advanced = false;
    char step_size_text[16] = "10000";
    bool step_size_edit = false;
    int refresh_counter = 0;
    float panel_scroll = 0.0f;  // Scroll offset for controls panel
    bool count_dist_edit = false;  // Dropdown state
    bool size_dist_edit = false;   // Dropdown state

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
        // Toggle maximize: Cmd+F on macOS, F11 on Windows/Linux
#ifdef __APPLE__
        bool toggle_maximize = (IsKeyDown(KEY_LEFT_SUPER) || IsKeyDown(KEY_RIGHT_SUPER)) && IsKeyPressed(KEY_F);
#else
        bool toggle_maximize = IsKeyPressed(KEY_F11);
#endif
        if (toggle_maximize) {
            if (IsWindowMaximized()) {
                RestoreWindow();
            } else {
                MaximizeWindow();
            }
        }

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

        // Update simulation
        if (running && !sim.stats.collapsed) {
            sim_run(&sim, (int)gens_per_frame);

            // Refresh colorizer cache periodically
            refresh_counter++;
            if (refresh_counter >= 60) {  // Every second at 60 FPS
                refresh_counter = 0;
                // Don't clear cache - let it grow for performance
            }
        }

        // Drawing
        BeginDrawing();
        ClearBackground((Color){30, 30, 35, 255});

        // Draw grid
        draw_grid(&sim, &colorizer, 10, 10, grid_width);

        // Draw panel background
        DrawRectangle(panel_x, 0, PANEL_WIDTH, screen_height, (Color){25, 25, 30, 255});

        // Handle scrolling when mouse is over panel
        Rectangle panel_rect = {panel_x, 0, PANEL_WIDTH, screen_height};
        Vector2 mouse = GetMousePosition();

        // Calculate content height (approximate based on controls)
        int content_height = 760;  // Base height (reduced since checkbox moved)
        if (show_advanced) {
            content_height += 170;
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
        DrawText("(Cmd+F toggle maximize)", panel_x + 20, 48, 12, GRAY);
#else
        DrawText("(F11 toggle maximize)", panel_x + 20, 48, 12, GRAY);
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
        }
        btn_y += btn_spacing;

        if (GuiButton((Rectangle){panel_x + 20, btn_y, 180, btn_h}, TextFormat("#79#Step %d", step_size))) {
            sim_run(&sim, step_size);
        }
        if (GuiButton((Rectangle){panel_x + 210, btn_y, 180, btn_h}, "#07#Export FASTA")) {
            // Use macOS native save dialog via osascript
            char cmd[512];
            snprintf(cmd, sizeof(cmd),
                "osascript -e 'POSIX path of (choose file name with prompt \"Save FASTA as:\" default name \"censim_gen%d.fasta\")' 2>/dev/null",
                sim.stats.generation);
            FILE *pipe = popen(cmd, "r");
            if (pipe) {
                char filepath[1024] = {0};
                if (fgets(filepath, sizeof(filepath), pipe)) {
                    // Remove trailing newline
                    filepath[strcspn(filepath, "\n")] = 0;
                    if (strlen(filepath) > 0) {
                        FILE *f = fopen(filepath, "w");
                        if (f) {
                            for (int i = 0; i < sim.array.num_units; i++) {
                                fprintf(f, ">repeat_%d\n%s\n", i + 1, sim.array.units[i]);
                            }
                            fclose(f);
                            printf("Exported %d repeats to %s\n", sim.array.num_units, filepath);
                        }
                    }
                }
                pclose(pipe);
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
            hover_text = "Expected INDELs per generation (Poisson lambda)";
            hover_rect = indel_row;
        }
        btn_y += row_h;

        // INDEL size
        Rectangle size_row = {panel_x, btn_y - 5, PANEL_WIDTH, row_h};
        DrawText("INDEL size:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%.1f", sim.params.indel_size_lambda),
            &sim.params.indel_size_lambda, 1.0f, 100.0f);
        if (CheckCollisionPointRec(mouse, size_row)) {
            hover_text = "Expected repeat units per INDEL (Poisson lambda)";
            hover_rect = size_row;
        }
        btn_y += row_h;

        // SNP rate
        Rectangle snp_row = {panel_x, btn_y - 5, PANEL_WIDTH, row_h};
        DrawText("SNP rate:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%.2f", sim.params.snp_rate),
            &sim.params.snp_rate, 0.0f, 1.0f);
        if (CheckCollisionPointRec(mouse, snp_row)) {
            hover_text = "Expected SNPs per generation (Poisson lambda)";
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
        if (GuiButton((Rectangle){panel_x + 20, btn_y, 370, btn_h}, "#27#UMAP Visualization (slow)")) {
            // Get output path via save dialog
            char cmd[512];
            snprintf(cmd, sizeof(cmd),
                "osascript -e 'POSIX path of (choose file name with prompt \"Save UMAP visualization as:\" default name \"umap_gen%d.png\")' 2>/dev/null",
                sim.stats.generation);
            FILE *pipe = popen(cmd, "r");
            if (pipe) {
                char outpath[1024] = {0};
                if (fgets(outpath, sizeof(outpath), pipe)) {
                    outpath[strcspn(outpath, "\n")] = 0;
                    if (strlen(outpath) > 0) {
                        // Write temp FASTA
                        char tempfasta[256];
                        snprintf(tempfasta, sizeof(tempfasta), "/tmp/censim_temp_%d.fasta", sim.stats.generation);
                        FILE *f = fopen(tempfasta, "w");
                        if (f) {
                            for (int i = 0; i < sim.array.num_units; i++) {
                                fprintf(f, ">repeat_%d\n%s\n", i + 1, sim.array.units[i]);
                            }
                            fclose(f);

                            // Build and run visualization command
                            char script_cmd[2048];
                            get_umap_command(script_cmd, sizeof(script_cmd), tempfasta, outpath, grid_width);
                            printf("Running: %s\n", script_cmd);
                            system(script_cmd);
                        }
                    }
                }
                pclose(pipe);
            }
        }
        btn_y += btn_spacing;

        // Stats panel
        draw_stats(&sim, panel_x + 10, btn_y, running);
        btn_y += 225;

        // Advanced options (collapsible)
        if (GuiButton((Rectangle){panel_x + 20, btn_y, 370, 25},
                      show_advanced ? "#120#Advanced Options" : "#119#Advanced Options")) {
            show_advanced = !show_advanced;
        }
        btn_y += 30;

        if (show_advanced) {
            // Calculate dynamic height based on visible parameter sliders
            int adv_height = 165;
            if (sim.params.count_dist == DIST_NEGATIVE_BINOMIAL) adv_height += 26;
            if (sim.params.size_dist == SIZE_POWER_LAW) adv_height += 26;

            DrawRectangle(panel_x + 10, btn_y, PANEL_WIDTH - 20, adv_height, (Color){40, 40, 40, 200});

            int adv_y = btn_y + 10;

            // Step size
            DrawText("Step size:", panel_x + 20, adv_y, 16, LIGHTGRAY);
            if (GuiTextBox((Rectangle){panel_x + 120, adv_y - 3, 100, 24},
                          step_size_text, 16, step_size_edit)) {
                step_size_edit = !step_size_edit;
            }
            if (!step_size_edit) {
                int val = atoi(step_size_text);
                if (val > 0) step_size = val;
            }
            adv_y += 32;

            // Dup/Del bias (own row)
            DrawText("Dup/Del bias:", panel_x + 20, adv_y, 16, LIGHTGRAY);
            GuiSlider((Rectangle){panel_x + 140, adv_y - 2, 200, 20}, NULL, NULL,
                      &sim.params.dup_bias, 0.0f, 1.0f);
            DrawText(TextFormat("%.0f%% dup", sim.params.dup_bias * 100), panel_x + 350, adv_y, 14, WHITE);
            adv_y += 32;

            // Hard bounds checkbox
            GuiCheckBox((Rectangle){panel_x + 20, adv_y, 20, 20}, "Hard bounds (min/max)", &sim.params.bounding_enabled);
            adv_y += 35;

            // Events dropdown row
            int events_y = adv_y;
            DrawText("Events:", panel_x + 20, events_y + 2, 16, LIGHTGRAY);
            adv_y += 28;

            // Dispersion slider row (only if NB selected)
            int dispersion_y = adv_y;
            if (sim.params.count_dist == DIST_NEGATIVE_BINOMIAL) {
                DrawText("dispersion:", panel_x + 40, dispersion_y + 2, 14, GRAY);
                GuiSlider((Rectangle){panel_x + 140, dispersion_y, 180, 18}, NULL, NULL,
                          &sim.params.nb_dispersion, 0.1f, 5.0f);
                DrawText(TextFormat("%.1f", sim.params.nb_dispersion), panel_x + 330, dispersion_y + 2, 12, WHITE);
                adv_y += 26;
            }

            // Sizes dropdown row
            int sizes_y = adv_y;
            DrawText("Sizes:", panel_x + 20, sizes_y + 2, 16, LIGHTGRAY);
            adv_y += 28;

            // Alpha slider row (only if power law selected)
            if (sim.params.size_dist == SIZE_POWER_LAW) {
                DrawText("alpha:", panel_x + 40, adv_y + 2, 14, GRAY);
                GuiSlider((Rectangle){panel_x + 140, adv_y, 180, 18}, NULL, NULL,
                          &sim.params.power_law_alpha, 1.5f, 4.0f);
                DrawText(TextFormat("%.1f", sim.params.power_law_alpha), panel_x + 330, adv_y + 2, 12, WHITE);
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

        EndDrawing();
    }

    // Cleanup
    sim_free(&sim);
    colorizer_free(&colorizer);
    CloseWindow();

    return 0;
}

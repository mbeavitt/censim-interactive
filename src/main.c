#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <raylib.h>

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

#include "config.h"
#include "simulation.h"
#include "colorizer.h"

// ============================================================================
// Grid rendering
// ============================================================================

static void draw_grid(Simulation *sim, Colorizer *colorizer, int offset_x, int offset_y) {
    int num_units = sim->array.num_units;
    if (num_units == 0) return;

    for (int i = 0; i < num_units; i++) {
        int row = i / GRID_WIDTH;
        int col = i % GRID_WIDTH;

        Color c = colorizer_get_color(colorizer, sim->array.units[i]);

        int x = offset_x + col * TILE_SIZE;
        int y = offset_y + row * TILE_SIZE;

        DrawRectangle(x, y, TILE_SIZE - 1, TILE_SIZE - 1, c);
    }
}

// ============================================================================
// Stats panel
// ============================================================================

static void draw_stats(Simulation *sim, int x, int y) {
    int unique = sim_count_unique(sim);
    float diversity = (sim->array.num_units > 0)
        ? (float)unique / (float)sim->array.num_units
        : 0.0f;

    char buf[256];

    DrawRectangle(x, y, PANEL_WIDTH - 20, 200, (Color){40, 40, 40, 220});
    DrawRectangleLines(x, y, PANEL_WIDTH - 20, 200, LIGHTGRAY);

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

    if (sim->stats.collapsed) {
        DrawText("COLLAPSED!", x + 10, line, 18, RED);
    }
}

// ============================================================================
// Main
// ============================================================================

int main(void) {
    // Initialize window
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Centromere Evolution Simulator");
    SetTargetFPS(60);

    // Start in fullscreen
    ToggleFullscreen();

    // Initialize simulation
    Simulation sim;
    sim_init(&sim, DEFAULT_INITIAL_SIZE);

    // Initialize colorizer
    Colorizer colorizer;
    colorizer_init(&colorizer, 42);

    // UI state
    bool running = false;
    float gens_per_frame = 100.0f;
    int refresh_counter = 0;

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
        // Toggle fullscreen with F11
        if (IsKeyPressed(KEY_F11) || IsKeyPressed(KEY_F)) {
            ToggleFullscreen();
        }

        // Get current screen dimensions
        int screen_width = GetScreenWidth();
        int screen_height = GetScreenHeight();
        int panel_x = screen_width - PANEL_WIDTH;

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
        draw_grid(&sim, &colorizer, 10, 10);

        // Draw panel background
        DrawRectangle(panel_x, 0, PANEL_WIDTH, screen_height, (Color){25, 25, 30, 255});

        // Title
        DrawText("Controls", panel_x + 20, 20, 24, WHITE);
        DrawText("(F11 toggle fullscreen)", panel_x + 20, 48, 12, GRAY);

        // Buttons
        int btn_y = 75;
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

        if (GuiButton((Rectangle){panel_x + 20, btn_y, 180, btn_h}, "#79#Step 1000")) {
            sim_run(&sim, 1000);
        }
        btn_y += btn_spacing + 20;

        // Sliders
        DrawText("Parameters", panel_x + 20, btn_y, 20, WHITE);
        btn_y += 30;

        int slider_w = 260;
        int slider_h = 20;
        int label_w = 100;

        // INDEL rate
        DrawText("INDEL rate:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%.2f", sim.params.indel_rate),
            &sim.params.indel_rate, 0.0f, 3.0f);
        btn_y += 35;

        // INDEL size
        DrawText("INDEL size:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%.1f", sim.params.indel_size_lambda),
            &sim.params.indel_size_lambda, 1.0f, 30.0f);
        btn_y += 35;

        // SNP rate
        DrawText("SNP rate:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%.2f", sim.params.snp_rate),
            &sim.params.snp_rate, 0.0f, 1.0f);
        btn_y += 35;

        // Gens per frame
        DrawText("Gens/frame:", panel_x + 20, btn_y + 2, 16, LIGHTGRAY);
        GuiSlider(
            (Rectangle){panel_x + label_w + 20, btn_y, slider_w, slider_h},
            NULL, TextFormat("%d", (int)gens_per_frame),
            &gens_per_frame, 10.0f, 100000.0f);
        btn_y += 45;

        // Bounding checkbox
        GuiCheckBox((Rectangle){panel_x + 20, btn_y, 20, 20}, "Bounding enabled", &sim.params.bounding_enabled);
        btn_y += 50;

        // Stats panel
        draw_stats(&sim, panel_x + 10, btn_y);

        // FPS counter
        DrawFPS(screen_width - 100, 10);

        EndDrawing();
    }

    // Cleanup
    sim_free(&sim);
    colorizer_free(&colorizer);
    CloseWindow();

    return 0;
}

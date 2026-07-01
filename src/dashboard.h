#ifndef DASHBOARD_H
#define DASHBOARD_H

#include <stdbool.h>
#include "batch.h"
#include "reference.h"

// Multi-trajectory dashboard view: configure + launch a batch of threaded
// trajectories and watch the report's metric distributions build up live as
// survivors land. Owns its Batch. Drawn as a mode toggled from the single view.

typedef struct {
    // Editable config (slider-backed floats mirror integer params).
    float f_num_traj;
    float f_initial;
    float f_target_gens;
    float f_collapse;
    float f_indel_rate;
    float f_snp_rate;
    float f_indel_size;
    float f_size_ratio;    // dup:del ratio: size split AND (inverse) frequency split, for zero net drift
    bool  elastic;         // false = free drift (default); true = elastic bounds toward start size
    float f_elasticity;    // elastic pull strength (when elastic enabled)

    float f_nbins;         // target on-screen bars per plot (adaptive; applied live)
    float f_cheb_order;    // target polynomial order for Chebyshev fit
    char fit_text[6][128]; // stores the text readout of the fit for each plot
    bool  plot_log_y[6];   // per-plot log-scale Y toggle
    bool  autoscale_x;     // fit each plot's x-window to the data (~99%)
    bool  show_advanced;   // advanced controls collapsible state

    Reference ref;         // reference/ghost distributions (may be unloaded)
    bool  show_ref;        // overlay the ghost at all
    bool  ghost_bars;      // draw the ghost as faded transparent bars
    bool  ghost_fit;       // draw the ghost's fitted curve (same fit as the batch)
    float ghost_alpha;     // ghost bar opacity, 0..1
    char  ghost_name[80];  // label for the current ghost source (file / pinned run)

    Batch batch;
    bool  has_batch;       // batch_init has been called (needs free)
    bool  started;         // a run is/has been launched

    // Parameter sweep state
    bool  sweep_enabled[6];
    float f_sweep_steps;
    float f_sweep_max_min;
    bool  sweep_running;
    int   sweep_total_runs;
    int   sweep_current_run;
    double sweep_batch_start_time;
    void *sweep_summary_file; // FILE*
    char  sweep_dir[256];
    int   sweep_indices[6];
    int   sweep_counts[6];
    bool  sweep_browsing;
    int   browse_run_idx;
    int   browse_run_count;   // number of run_*.csv in the browsed sweep (for the scrubber)
    float f_browse_run;       // slider-backed run index (1-based)
} Dashboard;

void dashboard_init(Dashboard *d);
void dashboard_free(Dashboard *d);

// Run one frame of the dashboard (controls + plots). `panel_w` is the width of
// the right-hand control panel; plots fill the area to its left.
void dashboard_update_draw(Dashboard *d, int screen_w, int screen_h, int panel_w);

#endif // DASHBOARD_H

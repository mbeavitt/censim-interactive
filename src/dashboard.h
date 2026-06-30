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

    Reference ref;         // real-data reference distributions (may be unloaded)
    bool  show_ref;        // overlay the real-data curve

    Batch batch;
    bool  has_batch;       // batch_init has been called (needs free)
    bool  started;         // a run is/has been launched
} Dashboard;

void dashboard_init(Dashboard *d);
void dashboard_free(Dashboard *d);

// Run one frame of the dashboard (controls + plots). `panel_w` is the width of
// the right-hand control panel; plots fill the area to its left.
void dashboard_update_draw(Dashboard *d, int screen_w, int screen_h, int panel_w);

#endif // DASHBOARD_H

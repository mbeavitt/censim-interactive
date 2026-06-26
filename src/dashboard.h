#ifndef DASHBOARD_H
#define DASHBOARD_H

#include <stdbool.h>
#include "batch.h"

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
    bool  unbounded;       // true = paper free-drift; false = elastic bounded

    float f_nbins;         // histogram bin resolution (applied on next Run)
    bool  plot_log_y[6];   // per-plot log-scale Y toggle
    bool  show_advanced;   // advanced controls collapsible state

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

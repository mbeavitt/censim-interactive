#ifndef STAINED_GLASS_H
#define STAINED_GLASS_H

// Live "stained glass" self-identity panel for the single-run view.
//
// Unit-vs-unit identity is approximated from MinHash sketches (Jaccard -> Mash
// identity via a precomputed lookup), so the O(N^2) matrix stays cheap. The
// array is subsampled to at most SG_MAX_UNITS units per refresh and binned into
// an SG_DIM x SG_DIM texture.
//
// The ~20 ms recompute runs on a background worker thread so it never stalls the
// 60 fps view: the main thread snapshots the array, hands it to the worker, and
// uploads the finished image (GPU calls stay main-thread) whenever it is ready.
// sg_update_draw throttles submissions to a target refresh interval.

#include <stdint.h>
#include <stdbool.h>
#include <pthread.h>

#include <raylib.h>
#include "simulation.h"

#define SG_KMER_K      16     // 16 bases pack into a uint32 (2 bits each)
#define SG_SKETCH      64     // MinHash sketch size
#define SG_DIM         384    // texture / matrix resolution (square)
#define SG_MAX_UNITS   2000   // units compared per refresh (subsampled from array)
#define SG_SEQ_STRIDE  256    // max stored length per snapshotted unit (monomer ~178)

// Producer/consumer states (guarded by mtx).
enum { SG_IDLE = 0, SG_REQUESTED, SG_COMPUTING, SG_DONE };

typedef struct {
    // --- main-thread only ---
    Texture2D tex;
    bool      computed;       // has a first result been uploaded
    double    last_submit;    // GetTime() of last work submission
    int       last_gen;       // sim generation at last submission
    int       last_units;     // array size at last submission
    float     lo, hi;         // last uploaded colour range (for the label)

    // --- shared, guarded by mtx ---
    pthread_t       thread;
    pthread_mutex_t mtx;
    pthread_cond_t  cv;
    int             state;
    bool            quit;
    char           *snap;     // SG_MAX_UNITS*SG_SEQ_STRIDE: main writes, worker reads
    int             snap_n;   // units in the current snapshot
    Color          *back;     // SG_DIM*SG_DIM: worker writes, main uploads
    float           res_lo, res_hi;

    // --- worker-private scratch ---
    uint32_t *sig;            // SG_MAX_UNITS*SG_SKETCH
    float    *avg;            // SG_DIM*SG_DIM
    uint32_t *cnt;            // SG_DIM*SG_DIM
    int      *bin;            // SG_MAX_UNITS
    float     ema_lo, ema_hi; // smoothed colour range (persists across refreshes)
    bool      ema_init;
} StainedGlass;

void sg_init(StainedGlass *sg);
void sg_free(StainedGlass *sg);

// Upload any finished result, then (if at least interval_s has elapsed and the
// array changed) snapshot the current array and kick a background recompute.
// Finally draw the panel into `box`. Call every frame from the main thread.
void sg_update_draw(StainedGlass *sg, const Simulation *sim, Rectangle box, double interval_s);

#endif // STAINED_GLASS_H

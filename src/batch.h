#ifndef BATCH_H
#define BATCH_H

#include <pthread.h>
#include <stdbool.h>
#include "simulation.h"
#include "hor.h"
#include "hist.h"

// Multi-trajectory batch runner.
//
// A pool of worker threads (= core count) runs `num_trajectories` independent
// simulations to `target_generations` (or until collapse). Each worker holds one
// array at a time, computes that trajectory's metrics, records a result, frees the
// array, and takes the next trajectory — so memory stays bounded at ~num_workers
// arrays regardless of batch size.
//
// Per the agreed data model: metric DISTRIBUTIONS are built from SURVIVORS only;
// collapsed trajectories contribute to the collapse rate and the time-to-collapse
// distribution. The UI thread reads the shared state under `lock` each frame.

typedef struct {
    int          num_trajectories;
    int          initial_size;
    long         target_generations;
    SimParams    base_params;     // mutation rates / distributions / bounding for every trajectory
    unsigned int seed_base;       // trajectory i uses seed_base + i
    int          nbins;           // histogram bin resolution (<=0 -> default 50)
} BatchConfig;

typedef struct {
    int    index;
    bool   done;
    bool   collapsed;
    long   generations;   // generations reached (target if survived, collapse gen if collapsed)
    int    final_units;
    // metric values (meaningful only if !collapsed)
    int    unique_units;
    double unique_per_kb;
    HorStats hor;
} TrajResult;

typedef struct {
    BatchConfig cfg;
    TrajResult *results;          // [num_trajectories]

    // Survivor distributions (one value per surviving trajectory)
    Histogram h_unique_per_kb;
    Histogram h_hors_per_kb;
    // Per-HOR distributions (many values per surviving trajectory)
    Histogram h_block_size;
    Histogram h_block_gap;
    Histogram h_similarity;
    Histogram h_diversity;
    Histogram h_composite;
    // Collapse timing (one value per collapsed trajectory)
    Histogram h_collapse_gen;

    // Progress (read by UI under lock)
    int  completed;
    int  survived;
    int  collapsed_count;

    // Internal scheduling / lifecycle
    int  next_index;              // work-queue cursor
    volatile bool stop_requested;
    bool running;
    int  num_workers;
    int  workers_running;         // live worker threads (for non-blocking stop/reap)
    pthread_t *threads;
    pthread_mutex_t lock;
} Batch;

// Detected hardware concurrency (>=1).
int batch_hardware_threads(void);

// Initialise a batch (allocates results + histograms, inits mutex). Does not start
// workers. `num_workers` <= 0 means "use batch_hardware_threads()".
void batch_init(Batch *b, BatchConfig cfg, int num_workers);

// Spawn the worker pool and begin running (non-blocking).
void batch_start(Batch *b);

// True once every trajectory has finished.
bool batch_is_complete(Batch *b);

// Non-blocking: signal workers to abandon all in-flight work and exit ASAP.
// The UI keeps running; reap with batch_join once workers_running hits 0.
void batch_request_stop(Batch *b);

// Request cancellation and join all workers (blocks until they exit).
void batch_stop(Batch *b);

// Join workers if the batch finished on its own (no-op if still running).
void batch_join(Batch *b);

// Free all resources. Stops workers first if still running.
void batch_free(Batch *b);

// Convenience: fill a SimParams with the paper baseline (unbounded drift), then the
// caller may tweak. collapse_threshold defaults to DEFAULT_COLLAPSE_THRESHOLD.
void batch_default_params(SimParams *p, bool unbounded);

#endif // BATCH_H

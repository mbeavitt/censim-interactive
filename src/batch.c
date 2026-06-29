#include "batch.h"
#include "config.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// Run the simulation in chunks so the worker can react to cancellation and check
// for collapse between chunks (sim_run itself stops early on collapse).
#define BATCH_CHUNK_GENS 50000

int batch_hardware_threads(void) {
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return (n > 0) ? (int)n : 1;
}

void batch_default_params(SimParams *p, bool unbounded) {
    p->indel_rate = DEFAULT_INDEL_RATE;          // Poisson 0.5
    p->indel_size_lambda = DEFAULT_INDEL_SIZE_LAMBDA;  // Poisson 7.6
    p->dup_del_size_ratio = DEFAULT_DUP_DEL_SIZE_RATIO; // 1.0 = dup/del equal size
    p->snp_rate = DEFAULT_SNP_RATE;              // Poisson 0.1
    p->min_array_size = DEFAULT_MIN_ARRAY_SIZE;
    p->max_array_size = DEFAULT_MAX_ARRAY_SIZE;
    p->collapse_threshold = DEFAULT_COLLAPSE_THRESHOLD;
    p->dup_bias = 0.5f;
    p->count_dist = DIST_POISSON;
    p->size_dist = SIZE_POISSON;
    p->nb_dispersion = 1.0f;
    p->power_law_alpha = 2.5f;
    if (unbounded) {
        // Paper regime: free drift, collapse when small.
        p->bounding_enabled = false;
        p->elasticity = 0.0f;
        p->target_size = DEFAULT_TRAJECTORY_INITIAL_SIZE;
    } else {
        // Bounded: soft elastic pull toward the starting size.
        p->bounding_enabled = false;
        p->elasticity = 0.15f;
        p->target_size = DEFAULT_TRAJECTORY_INITIAL_SIZE;
    }
}

void batch_init(Batch *b, BatchConfig cfg, int num_workers) {
    memset(b, 0, sizeof(*b));
    b->cfg = cfg;
    b->results = (TrajResult *)calloc(cfg.num_trajectories, sizeof(TrajResult));
    for (int i = 0; i < cfg.num_trajectories; i++) b->results[i].index = i;

    // Fixed binning ranges. The unbounded metrics use generous LOG ceilings so
    // even extreme-parameter runs stay in-range (a high ceiling is cheap on a log
    // axis, and the dashboard autoscales the display to the populated region).
    // These cover anything actually computable -- arrays large enough to exceed
    // them would make the O(n^2) HOR scan intractable anyway. unique/kb and
    // similarity have true upper bounds (~5.6 and 1) so stay linear/tight.
    int nb = (cfg.nbins > 0) ? cfg.nbins : 50;
    hist_init(&b->h_unique_per_kb, 0.0f, 6.0f, nb, 0);      // linear; hard max ~5.6 (all-unique)
    hist_init(&b->h_hors_per_kb,   1.0f, 1.0e6f, nb, 1);    // log X
    hist_init(&b->h_block_size,    1.0f, 1.0e6f, nb, 1);    // log X (block <= array size)
    hist_init(&b->h_block_gap,     1.0f, 1.0e7f, nb, 1);    // log X
    hist_init(&b->h_similarity,    0.0f, 1.0f, nb, 0);      // linear; bounded 0..1
    hist_init(&b->h_diversity,     0.0f, 1.0f, nb, 0);      // linear; bounded 0..1
    hist_init(&b->h_composite,     1.0f, 1.0e15f, nb, 1);   // log X
    float gmax = (cfg.target_generations > 0) ? (float)cfg.target_generations : 1.0f;
    hist_init(&b->h_collapse_gen,  0.0f, gmax, nb, 0);

    b->next_index = 0;
    b->stop_requested = false;
    b->running = false;
    b->num_workers = (num_workers > 0) ? num_workers : batch_hardware_threads();
    if (b->num_workers > cfg.num_trajectories && cfg.num_trajectories > 0)
        b->num_workers = cfg.num_trajectories;
    b->threads = (pthread_t *)calloc(b->num_workers, sizeof(pthread_t));
    pthread_mutex_init(&b->lock, NULL);
}

// Claim the next trajectory index, or -1 if none left / cancelled.
static int claim_index(Batch *b) {
    int idx = -1;
    pthread_mutex_lock(&b->lock);
    if (!b->stop_requested && b->next_index < b->cfg.num_trajectories)
        idx = b->next_index++;
    pthread_mutex_unlock(&b->lock);
    return idx;
}

static void *worker_main(void *arg) {
    Batch *b = (Batch *)arg;
    int idx;
    while ((idx = claim_index(b)) >= 0) {
        Simulation sim;
        sim_init(&sim, b->cfg.initial_size, b->cfg.seed_base + (unsigned int)idx);
        sim.params = b->cfg.base_params;  // apply batch template; rng_state stays as seeded by sim_init

        long target = b->cfg.target_generations;
        while (sim.stats.generation < target && !sim.stats.collapsed && !b->stop_requested) {
            long remaining = target - sim.stats.generation;
            int chunk = (remaining < BATCH_CHUNK_GENS) ? (int)remaining : BATCH_CHUNK_GENS;
            sim_run(&sim, chunk);
        }

        bool reached  = (sim.stats.generation >= target);
        bool collapsed = sim.stats.collapsed;
        bool survived  = reached && !collapsed;

        TrajResult r;
        memset(&r, 0, sizeof(r));
        r.index = idx;
        r.done = true;
        r.collapsed = collapsed;
        r.generations = sim.stats.generation;
        r.final_units = sim.array.num_units;

        if (survived) {
            r.unique_units = sim_count_unique(&sim);
            double kb = (double)sim.array.num_units * REPEAT_SIZE / 1000.0;
            r.unique_per_kb = (kb > 0.0) ? r.unique_units / kb : 0.0;

            // Fold per-HOR values into thread-local histograms (no lock during the
            // O(n^2) scan), then merge into the shared ones under lock.
            Histogram ls, lg, lsi, ld, lc;
            hist_init(&ls,  b->h_block_size.min, b->h_block_size.max, b->h_block_size.nbins, b->h_block_size.log_scale);
            hist_init(&lg,  b->h_block_gap.min,  b->h_block_gap.max,  b->h_block_gap.nbins,  b->h_block_gap.log_scale);
            hist_init(&lsi, b->h_similarity.min, b->h_similarity.max, b->h_similarity.nbins, b->h_similarity.log_scale);
            hist_init(&ld,  b->h_diversity.min,  b->h_diversity.max,  b->h_diversity.nbins,  b->h_diversity.log_scale);
            hist_init(&lc,  b->h_composite.min,  b->h_composite.max,  b->h_composite.nbins,  b->h_composite.log_scale);
            HorHistSet hs = { &ls, &lg, &lsi, &ld, &lc };
            hor_scan(&sim.array, &r.hor, &hs, &b->stop_requested);

            // If stop was requested, the scan may have aborted mid-way; discard.
            if (!b->stop_requested) {
                pthread_mutex_lock(&b->lock);
                b->results[idx] = r;
                hist_add(&b->h_unique_per_kb, (float)r.unique_per_kb);
                hist_add(&b->h_hors_per_kb,   (float)r.hor.hors_per_kb);
                hist_merge(&b->h_block_size, &ls);
                hist_merge(&b->h_block_gap,  &lg);
                hist_merge(&b->h_similarity, &lsi);
                hist_merge(&b->h_diversity,  &ld);
                hist_merge(&b->h_composite,  &lc);
                b->survived++;
                b->completed++;
                pthread_mutex_unlock(&b->lock);
            }

            hist_free(&ls); hist_free(&lg); hist_free(&lsi); hist_free(&ld); hist_free(&lc);
        } else {
            pthread_mutex_lock(&b->lock);
            b->results[idx] = r;
            if (collapsed) {
                hist_add(&b->h_collapse_gen, (float)r.generations);
                b->collapsed_count++;
                b->completed++;
            }
            // aborted-by-stop trajectories are left uncounted in completed
            pthread_mutex_unlock(&b->lock);
        }

        sim_free(&sim);
    }
    pthread_mutex_lock(&b->lock);
    b->workers_running--;
    pthread_mutex_unlock(&b->lock);
    return NULL;
}

void batch_start(Batch *b) {
    b->running = true;
    b->workers_running = b->num_workers;
    for (int i = 0; i < b->num_workers; i++) {
        pthread_create(&b->threads[i], NULL, worker_main, b);
    }
}

void batch_request_stop(Batch *b) {
    b->stop_requested = true;  // non-blocking; workers abandon work and exit
}

bool batch_is_complete(Batch *b) {
    pthread_mutex_lock(&b->lock);
    bool done = (b->completed >= b->cfg.num_trajectories);
    pthread_mutex_unlock(&b->lock);
    return done;
}

void batch_join(Batch *b) {
    if (!b->running) return;
    for (int i = 0; i < b->num_workers; i++) pthread_join(b->threads[i], NULL);
    b->running = false;
}

void batch_stop(Batch *b) {
    b->stop_requested = true;
    batch_join(b);
}

void batch_free(Batch *b) {
    if (b->running) batch_stop(b);
    hist_free(&b->h_unique_per_kb);
    hist_free(&b->h_hors_per_kb);
    hist_free(&b->h_block_size);
    hist_free(&b->h_block_gap);
    hist_free(&b->h_similarity);
    hist_free(&b->h_diversity);
    hist_free(&b->h_composite);
    hist_free(&b->h_collapse_gen);
    free(b->results);
    free(b->threads);
    pthread_mutex_destroy(&b->lock);
    memset(b, 0, sizeof(*b));
}

// Single-threaded benchmark runner for the simulation + HOR detector.
//
// Deliberately NOT the GUI app and NOT the threaded batch runner: it runs
// trajectories serially so a sampling profiler (perf / samply / gperftools /
// callgrind) sees a clean call tree of the actual hot paths (the O(n^2) HOR scan
// and the sim mutation loop), with no raylib or pthread noise.
//
// Build:   premake5 gmake2 && make config=profile
// Run:     ./bin/Profile/bench [options]
// Profile: CPUPROFILE=cpu.prof LD_PRELOAD=libprofiler.so ./bin/Profile/bench --scan-only 20
//
// Options (all optional):
//   --traj N        number of trajectories to run             (default 8)
//   --gens G        generations per trajectory                (default 1000000)
//   --size S        starting array size in units              (default 1500)
//   --collapse C    collapse threshold in units               (default 300)
//   --bounded       elastic-bounded instead of free drift     (default unbounded)
//   --seed K        base RNG seed                             (default 1)
//   --scan-only R   build ONE array (to --gens) then run the HOR scan R times,
//                   isolating the detector from the simulation
//   --no-scan       run the simulation only (skip HOR scanning)
#include "simulation.h"
#include "hor.h"
#include "hist.h"
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double cpu_s(void) { return (double)clock() / (double)CLOCKS_PER_SEC; }

static void set_params(SimParams *p, int unbounded, int collapse, int size) {
    p->indel_rate        = DEFAULT_INDEL_RATE;          // Poisson 0.5
    p->indel_size_lambda = DEFAULT_INDEL_SIZE_LAMBDA;   // Poisson 7.6
    p->snp_rate          = DEFAULT_SNP_RATE;            // Poisson 0.1
    p->min_array_size    = DEFAULT_MIN_ARRAY_SIZE;
    p->max_array_size    = DEFAULT_MAX_ARRAY_SIZE;
    p->collapse_threshold = collapse;
    p->dup_bias          = 0.5f;
    p->count_dist        = DIST_POISSON;
    p->size_dist         = SIZE_POISSON;
    p->nb_dispersion     = 1.0f;
    p->power_law_alpha   = 2.5f;
    p->bounding_enabled  = false;
    p->elasticity        = unbounded ? 0.0f : 0.15f;
    p->target_size       = size;
}

static int arg_int(int argc, char **argv, const char *key, int def) {
    for (int i = 1; i < argc - 1; i++)
        if (strcmp(argv[i], key) == 0) return atoi(argv[i + 1]);
    return def;
}
static int arg_flag(int argc, char **argv, const char *key) {
    for (int i = 1; i < argc; i++) if (strcmp(argv[i], key) == 0) return 1;
    return 0;
}

int main(int argc, char **argv) {
    int traj      = arg_int(argc, argv, "--traj", 8);
    long gens     = arg_int(argc, argv, "--gens", 1000000);
    int size      = arg_int(argc, argv, "--size", 1500);
    int collapse  = arg_int(argc, argv, "--collapse", 300);
    unsigned seed = (unsigned)arg_int(argc, argv, "--seed", 1);
    int unbounded = !arg_flag(argc, argv, "--bounded");
    int scan_only = arg_int(argc, argv, "--scan-only", 0);
    int no_scan   = arg_flag(argc, argv, "--no-scan");

    printf("bench: traj=%d gens=%ld size=%d collapse=%d %s%s%s\n",
           traj, gens, size, collapse, unbounded ? "unbounded" : "bounded",
           no_scan ? " no-scan" : "", scan_only ? " scan-only" : "");

    // --- scan-only: isolate the HOR detector on a single matured array ---
    if (scan_only > 0) {
        Simulation sim;
        sim_init(&sim, size, seed);
        set_params(&sim.params, unbounded, collapse, size);
        sim_run(&sim, gens);
        printf("array: %d units after %ld gens (collapsed=%d)\n",
               sim.array.num_units, gens, sim.stats.collapsed);
        double t0 = cpu_s();
        long total_hors = 0;
        for (int r = 0; r < scan_only; r++) {
            HorStats s;
            hor_scan(&sim.array, &s, NULL, NULL);
            total_hors = s.num_hors;
        }
        double dt = cpu_s() - t0;
        printf("HOR scan x%d: %.3fs total, %.4fs/scan  (%ld HORs, %.1f Mpairs/s)\n",
               scan_only, dt, dt / scan_only, total_hors,
               (double)sim.array.num_units * sim.array.num_units / 2.0 / (dt / scan_only) / 1e6);
        sim_free(&sim);
        return 0;
    }

    // --- full serial sweep: sim + (optional) scan over `traj` trajectories ---
    double sim_t = 0, scan_t = 0;
    int survived = 0, collapsed = 0;
    double sum_uniq_per_kb = 0, sum_hors_per_kb = 0;

    for (int i = 0; i < traj; i++) {
        Simulation sim;
        sim_init(&sim, size, seed + (unsigned)i);
        set_params(&sim.params, unbounded, collapse, size);

        double t0 = cpu_s();
        sim_run(&sim, gens);
        sim_t += cpu_s() - t0;

        int reached = (sim.stats.generation >= gens);
        if (reached && !sim.stats.collapsed) {
            survived++;
            if (!no_scan) {
                HorStats s;
                double t1 = cpu_s();
                hor_scan(&sim.array, &s, NULL, NULL);
                scan_t += cpu_s() - t1;
                double kb = (double)sim.array.num_units * REPEAT_SIZE / 1000.0;
                sum_uniq_per_kb += kb > 0 ? sim_count_unique(&sim) / kb : 0;
                sum_hors_per_kb += s.hors_per_kb;
            }
        } else if (sim.stats.collapsed) {
            collapsed++;
        }
        sim_free(&sim);
    }

    printf("survived=%d collapsed=%d\n", survived, collapsed);
    if (survived && !no_scan) {
        printf("means (survivors): unique/kb=%.2f  HORs/kb=%.1f\n",
               sum_uniq_per_kb / survived, sum_hors_per_kb / survived);
    }
    printf("time: sim=%.3fs  scan=%.3fs  total=%.3fs\n", sim_t, scan_t, sim_t + scan_t);
    return 0;
}

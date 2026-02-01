#ifndef SIMULATION_H
#define SIMULATION_H

#include <stdbool.h>

// Distribution types for mutation modeling
typedef enum {
    DIST_POISSON = 0,       // Mean = variance, standard choice
    DIST_NEGATIVE_BINOMIAL, // Overdispersed (variance > mean), common in biology
    DIST_COUNT              // Number of distribution types
} CountDistribution;

typedef enum {
    SIZE_POISSON = 0,   // Symmetric around mean
    SIZE_GEOMETRIC,     // Exponential decay (small events more likely)
    SIZE_POWER_LAW,     // Heavy tail (rare large events)
    SIZE_COUNT          // Number of size distribution types
} SizeDistribution;

// Forward declarations
typedef struct RepeatArray RepeatArray;
typedef struct SimParams SimParams;
typedef struct SimStats SimStats;
typedef struct Simulation Simulation;

// Repeat array - dynamic array of 178-char sequences
struct RepeatArray {
    char **units;       // Array of pointers to sequences
    int num_units;      // Current count
    int capacity;       // Allocated capacity
};

// Simulation parameters
struct SimParams {
    float indel_rate;
    float indel_size_lambda;
    float snp_rate;
    int min_array_size;
    int max_array_size;
    bool bounding_enabled;
    // Elastic bounding
    int target_size;
    float elasticity;  // 0 = no effect, higher = stronger pull toward target
    float dup_bias;    // 0 = all deletions, 0.5 = equal, 1 = all duplications
    // Distribution models
    CountDistribution count_dist;  // For event counts (SNPs, indels)
    SizeDistribution size_dist;    // For indel sizes
    float nb_dispersion;           // Dispersion parameter for negative binomial (higher = more overdispersion)
    float power_law_alpha;         // Shape parameter for power law (lower = heavier tail)
};

// Statistics
struct SimStats {
    int generation;
    int snp_count;
    int dup_count;
    int del_count;
    bool collapsed;
};

// Main simulation struct
struct Simulation {
    RepeatArray array;
    SimParams params;
    SimStats stats;
};

// Public API
void sim_init(Simulation *sim, int initial_size);
void sim_free(Simulation *sim);
void sim_step(Simulation *sim);
void sim_run(Simulation *sim, int generations);
void sim_reset(Simulation *sim);

// Statistics
int sim_count_unique(Simulation *sim);
float sim_diversity(Simulation *sim);

#endif // SIMULATION_H

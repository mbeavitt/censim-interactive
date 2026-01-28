#ifndef SIMULATION_H
#define SIMULATION_H

#include <stdbool.h>

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

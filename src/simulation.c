#include "simulation.h"
#include "config.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// ============================================================================
// Random number generation
// ============================================================================

static unsigned int rng_state = 42;

static void rng_seed(unsigned int seed) {
    rng_state = seed;
}

// Xorshift32
static unsigned int rng_next(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

static float rng_float(void) {
    return (float)(rng_next() & 0x7FFFFFFF) / (float)0x7FFFFFFF;
}

static int rng_int(int max) {
    return (int)(rng_next() % (unsigned int)max);
}

// Poisson sampling using inverse transform
static int sample_poisson(float lambda) {
    if (lambda <= 0) return 0;

    float L = expf(-lambda);
    int k = 0;
    float p = 1.0f;

    do {
        k++;
        p *= rng_float();
    } while (p > L);

    return k - 1;
}

// ============================================================================
// Repeat Array operations
// ============================================================================

static void array_init(RepeatArray *arr, int capacity) {
    arr->capacity = capacity;
    arr->num_units = 0;
    arr->units = (char **)malloc(capacity * sizeof(char *));
}

static void array_free(RepeatArray *arr) {
    for (int i = 0; i < arr->num_units; i++) {
        free(arr->units[i]);
    }
    free(arr->units);
    arr->units = NULL;
    arr->num_units = 0;
    arr->capacity = 0;
}

static void array_ensure_capacity(RepeatArray *arr, int needed) {
    if (needed <= arr->capacity) return;

    int new_capacity = arr->capacity;
    while (new_capacity < needed) {
        new_capacity *= 2;
    }

    arr->units = (char **)realloc(arr->units, new_capacity * sizeof(char *));
    arr->capacity = new_capacity;
}

static char *alloc_unit(const char *src) {
    char *unit = (char *)malloc(REPEAT_SIZE + 1);
    memcpy(unit, src, REPEAT_SIZE);
    unit[REPEAT_SIZE] = '\0';
    return unit;
}

// Tandem duplication: duplicate units[start:end] in place after end
static void duplicate_units(RepeatArray *arr, int start, int end) {
    int count = end - start;
    array_ensure_capacity(arr, arr->num_units + count);

    // Shift everything after 'end' to make room
    memmove(&arr->units[end + count], &arr->units[end],
            (arr->num_units - end) * sizeof(char *));

    // Copy the units
    for (int i = 0; i < count; i++) {
        arr->units[end + i] = alloc_unit(arr->units[start + i]);
    }

    arr->num_units += count;
}

// Delete units[start:end]
static void delete_units(RepeatArray *arr, int start, int end) {
    int count = end - start;

    // Free the deleted units
    for (int i = start; i < end; i++) {
        free(arr->units[i]);
    }

    // Shift remaining units
    memmove(&arr->units[start], &arr->units[end],
            (arr->num_units - end) * sizeof(char *));

    arr->num_units -= count;
}

// ============================================================================
// Mutation functions
// ============================================================================

static void apply_snps(Simulation *sim) {
    int n_snps = sample_poisson(sim->params.snp_rate);

    for (int i = 0; i < n_snps; i++) {
        if (sim->array.num_units == 0) break;

        // Pick random unit and position
        int unit_idx = rng_int(sim->array.num_units);
        int pos = rng_int(REPEAT_SIZE);

        char *unit = sim->array.units[unit_idx];
        char old_base = unit[pos];

        // Pick a different base
        char new_base;
        do {
            new_base = BASES[rng_int(4)];
        } while (new_base == old_base);

        unit[pos] = new_base;
        sim->stats.snp_count++;
    }
}

static void apply_indels(Simulation *sim) {
    int n_indels = sample_poisson(sim->params.indel_rate);

    for (int i = 0; i < n_indels; i++) {
        if (sim->array.num_units == 0) {
            sim->stats.collapsed = true;
            return;
        }

        // Elastic bounding: bias dup/del probability based on distance from target
        float dup_prob = 0.5f;
        if (sim->params.elasticity > 0.0f) {
            // How far are we from target? Positive = too big, negative = too small
            float deviation = (float)(sim->array.num_units - sim->params.target_size)
                            / (float)sim->params.target_size;
            // Bias: when too big, reduce dup probability; when too small, increase it
            // dup_prob = 0.5 - elasticity * deviation (clamped to [0.1, 0.9])
            dup_prob = 0.5f - sim->params.elasticity * deviation;
            if (dup_prob < 0.1f) dup_prob = 0.1f;
            if (dup_prob > 0.9f) dup_prob = 0.9f;
        }

        // Choose duplication or deletion based on biased probability
        bool is_dup = rng_float() < dup_prob;

        // Sample size
        int indel_size = sample_poisson(sim->params.indel_size_lambda);
        if (indel_size < 1) indel_size = 1;

        // Pick start position
        int start = rng_int(sim->array.num_units);
        int end = start + indel_size;

        // Bounds check
        if (end > sim->array.num_units) {
            continue;  // Skip this INDEL
        }

        if (is_dup) {
            // Check max size
            if (sim->params.bounding_enabled &&
                sim->array.num_units + indel_size > sim->params.max_array_size) {
                continue;
            }
            duplicate_units(&sim->array, start, end);
            sim->stats.dup_count++;
        } else {
            // Check min size
            if (sim->params.bounding_enabled &&
                sim->array.num_units - indel_size < sim->params.min_array_size) {
                continue;
            }
            delete_units(&sim->array, start, end);
            sim->stats.del_count++;
        }
    }

    // Check for collapse
    if (sim->array.num_units < sim->params.min_array_size) {
        sim->stats.collapsed = true;
    }
}

// ============================================================================
// Public API
// ============================================================================

void sim_init(Simulation *sim, int initial_size) {
    rng_seed((unsigned int)time(NULL));

    // Initialize array
    array_init(&sim->array, initial_size * 2);  // Extra capacity

    // Fill with copies of default monomer
    for (int i = 0; i < initial_size; i++) {
        sim->array.units[i] = alloc_unit(DEFAULT_MONOMER);
    }
    sim->array.num_units = initial_size;

    // Default parameters
    sim->params.indel_rate = DEFAULT_INDEL_RATE;
    sim->params.indel_size_lambda = DEFAULT_INDEL_SIZE_LAMBDA;
    sim->params.snp_rate = DEFAULT_SNP_RATE;
    sim->params.min_array_size = DEFAULT_MIN_ARRAY_SIZE;
    sim->params.max_array_size = DEFAULT_MAX_ARRAY_SIZE;
    sim->params.bounding_enabled = true;
    sim->params.target_size = DEFAULT_TARGET_SIZE;
    sim->params.elasticity = DEFAULT_ELASTICITY;

    // Reset stats
    sim->stats.generation = 0;
    sim->stats.snp_count = 0;
    sim->stats.dup_count = 0;
    sim->stats.del_count = 0;
    sim->stats.collapsed = false;
}

void sim_free(Simulation *sim) {
    array_free(&sim->array);
}

void sim_step(Simulation *sim) {
    if (sim->stats.collapsed) return;

    sim->stats.generation++;
    apply_snps(sim);
    apply_indels(sim);
}

void sim_run(Simulation *sim, int generations) {
    for (int i = 0; i < generations; i++) {
        if (sim->stats.collapsed) break;
        sim_step(sim);
    }
}

void sim_reset(Simulation *sim) {
    // Free existing array
    array_free(&sim->array);

    // Reinitialize
    int initial_size = DEFAULT_INITIAL_SIZE;
    array_init(&sim->array, initial_size * 2);

    for (int i = 0; i < initial_size; i++) {
        sim->array.units[i] = alloc_unit(DEFAULT_MONOMER);
    }
    sim->array.num_units = initial_size;

    // Reset stats but keep params
    sim->stats.generation = 0;
    sim->stats.snp_count = 0;
    sim->stats.dup_count = 0;
    sim->stats.del_count = 0;
    sim->stats.collapsed = false;
}

// Simple hash-based unique counting using FNV-1a
static unsigned int hash_sequence(const char *seq) {
    unsigned int hash = 2166136261u;
    for (int i = 0; i < REPEAT_SIZE; i++) {
        hash ^= (unsigned char)seq[i];
        hash *= 16777619u;
    }
    return hash;
}

int sim_count_unique(Simulation *sim) {
    if (sim->array.num_units == 0) return 0;

    // Use a simple hash set
    int table_size = sim->array.num_units * 2;
    unsigned int *seen = (unsigned int *)calloc(table_size, sizeof(unsigned int));
    int unique = 0;

    for (int i = 0; i < sim->array.num_units; i++) {
        unsigned int hash = hash_sequence(sim->array.units[i]);
        unsigned int idx = hash % table_size;

        // Linear probing
        while (seen[idx] != 0) {
            if (seen[idx] == hash) break;  // Collision assumed same (imperfect)
            idx = (idx + 1) % table_size;
        }

        if (seen[idx] == 0) {
            seen[idx] = hash;
            unique++;
        }
    }

    free(seen);
    return unique;
}

float sim_diversity(Simulation *sim) {
    if (sim->array.num_units == 0) return 0.0f;
    return (float)sim_count_unique(sim) / (float)sim->array.num_units;
}

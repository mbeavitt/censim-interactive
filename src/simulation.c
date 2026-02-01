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

// Gamma sampling using Marsaglia and Tsang's method
// Used internally for negative binomial
static float sample_gamma(float shape, float scale) {
    if (shape < 1.0f) {
        // For shape < 1, use shape+1 and then scale
        float u = rng_float();
        return sample_gamma(shape + 1.0f, scale) * powf(u, 1.0f / shape);
    }

    float d = shape - 1.0f / 3.0f;
    float c = 1.0f / sqrtf(9.0f * d);

    while (1) {
        float x, v;
        do {
            // Box-Muller for standard normal
            float u1 = rng_float();
            float u2 = rng_float();
            x = sqrtf(-2.0f * logf(u1 + 1e-10f)) * cosf(2.0f * 3.14159265f * u2);
            v = 1.0f + c * x;
        } while (v <= 0.0f);

        v = v * v * v;
        float u = rng_float();

        if (u < 1.0f - 0.0331f * x * x * x * x) {
            return d * v * scale;
        }

        if (logf(u + 1e-10f) < 0.5f * x * x + d * (1.0f - v + logf(v))) {
            return d * v * scale;
        }
    }
}

// Negative Binomial sampling
// Models overdispersed count data (variance > mean)
// Mean = mu, Variance = mu + mu^2/dispersion
// Higher dispersion = more overdispersion
static int sample_negative_binomial(float mu, float dispersion) {
    if (mu <= 0) return 0;
    if (dispersion <= 0) return sample_poisson(mu);  // Fall back to Poisson

    // NB as Gamma-Poisson mixture
    // r = dispersion, p = dispersion/(dispersion + mu)
    float r = dispersion;
    float lambda = sample_gamma(r, mu / r);
    return sample_poisson(lambda);
}

// Geometric sampling (number of failures before first success)
// Mean = (1-p)/p, so p = 1/(1+mean)
// Models exponentially decaying probability - small values more likely
static int sample_geometric(float mean) {
    if (mean <= 0) return 1;

    float p = 1.0f / (1.0f + mean);
    float u = rng_float();

    // Inverse transform: floor(log(u) / log(1-p))
    int result = (int)(logf(u + 1e-10f) / logf(1.0f - p + 1e-10f));
    return result < 1 ? 1 : result;
}

// Power law (Pareto) sampling
// P(X >= x) ~ x^(-alpha), heavy tailed
// Mean parameter controls the scale (x_min), alpha controls tail heaviness
// Small events common, but occasional large events occur
static int sample_power_law(float mean, float alpha) {
    if (mean <= 0) return 1;
    if (alpha <= 1.0f) alpha = 1.1f;  // Must be > 1 for finite mean

    // For Pareto: mean = alpha * x_min / (alpha - 1) when alpha > 1
    // So x_min = mean * (alpha - 1) / alpha
    float x_min = mean * (alpha - 1.0f) / alpha;
    if (x_min < 1.0f) x_min = 1.0f;

    float u = rng_float();

    // Inverse transform: x = x_min * (1 - u)^(-1/alpha)
    int result = (int)(x_min * powf(1.0f - u + 1e-10f, -1.0f / alpha));
    return result < 1 ? 1 : result;
}

// ============================================================================
// Distribution dispatch helpers
// ============================================================================

static int sample_count(CountDistribution dist, float mean, float dispersion) {
    switch (dist) {
        case DIST_NEGATIVE_BINOMIAL:
            return sample_negative_binomial(mean, dispersion);
        case DIST_POISSON:
        default:
            return sample_poisson(mean);
    }
}

static int sample_size(SizeDistribution dist, float mean, float power_law_alpha) {
    switch (dist) {
        case SIZE_GEOMETRIC:
            return sample_geometric(mean);
        case SIZE_POWER_LAW:
            return sample_power_law(mean, power_law_alpha);
        case SIZE_POISSON:
        default: {
            int s = sample_poisson(mean);
            return s < 1 ? 1 : s;
        }
    }
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
    int n_snps = sample_count(sim->params.count_dist, sim->params.snp_rate, sim->params.nb_dispersion);

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
    int n_indels = sample_count(sim->params.count_dist, sim->params.indel_rate, sim->params.nb_dispersion);

    for (int i = 0; i < n_indels; i++) {
        if (sim->array.num_units == 0) {
            sim->stats.collapsed = true;
            return;
        }

        // Start with base dup/del bias
        float dup_prob = sim->params.dup_bias;

        // Elastic bounding: further bias based on distance from target
        if (sim->params.elasticity > 0.0f) {
            // How far are we from target? Positive = too big, negative = too small
            float deviation = (float)(sim->array.num_units - sim->params.target_size)
                            / (float)sim->params.target_size;
            // Bias: when too big, reduce dup probability; when too small, increase it
            dup_prob = dup_prob - sim->params.elasticity * deviation;
            if (dup_prob < 0.05f) dup_prob = 0.05f;
            if (dup_prob > 0.95f) dup_prob = 0.95f;
        }

        // Choose duplication or deletion based on biased probability
        bool is_dup = rng_float() < dup_prob;

        // Sample size using selected distribution
        int indel_size = sample_size(sim->params.size_dist, sim->params.indel_size_lambda, sim->params.power_law_alpha);

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
    sim->params.dup_bias = 0.5f;  // Equal dup/del by default
    // Distribution defaults
    sim->params.count_dist = DIST_POISSON;
    sim->params.size_dist = SIZE_POISSON;
    sim->params.nb_dispersion = 1.0f;  // Moderate overdispersion when NB is used
    sim->params.power_law_alpha = 2.5f;  // Shape parameter (lower = heavier tail)

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

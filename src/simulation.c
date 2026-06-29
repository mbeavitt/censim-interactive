#include "simulation.h"
#include "config.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// ============================================================================
// Random number generation
// ============================================================================

// RNG state is per-trajectory (carried in Simulation, passed by pointer) so that
// trajectories run independently across threads with no shared mutable state.

// Xorshift32
static unsigned int rng_next(unsigned int *state) {
    unsigned int x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

static float rng_float(unsigned int *state) {
    return (float)(rng_next(state) & 0x7FFFFFFF) / (float)0x7FFFFFFF;
}

static int rng_int(unsigned int *state, int max) {
    return (int)(rng_next(state) % (unsigned int)max);
}

// Poisson sampling using inverse transform
static int sample_poisson(unsigned int *state, float lambda) {
    if (lambda <= 0) return 0;

    float L = expf(-lambda);
    int k = 0;
    float p = 1.0f;

    do {
        k++;
        p *= rng_float(state);
    } while (p > L);

    return k - 1;
}

// Gamma sampling using Marsaglia and Tsang's method
// Used internally for negative binomial
static float sample_gamma(unsigned int *state, float shape, float scale) {
    if (shape < 1.0f) {
        // For shape < 1, use shape+1 and then scale
        float u = rng_float(state);
        return sample_gamma(state, shape + 1.0f, scale) * powf(u, 1.0f / shape);
    }

    float d = shape - 1.0f / 3.0f;
    float c = 1.0f / sqrtf(9.0f * d);

    while (1) {
        float x, v;
        do {
            // Box-Muller for standard normal
            float u1 = rng_float(state);
            float u2 = rng_float(state);
            x = sqrtf(-2.0f * logf(u1 + 1e-10f)) * cosf(2.0f * 3.14159265f * u2);
            v = 1.0f + c * x;
        } while (v <= 0.0f);

        v = v * v * v;
        float u = rng_float(state);

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
static int sample_negative_binomial(unsigned int *state, float mu, float dispersion) {
    if (mu <= 0) return 0;
    if (dispersion <= 0) return sample_poisson(state, mu);  // Fall back to Poisson

    // NB as Gamma-Poisson mixture
    // r = dispersion, p = dispersion/(dispersion + mu)
    float r = dispersion;
    float lambda = sample_gamma(state, r, mu / r);
    return sample_poisson(state, lambda);
}

// Geometric sampling (number of failures before first success)
// Mean = (1-p)/p, so p = 1/(1+mean)
// Models exponentially decaying probability - small values more likely
static int sample_geometric(unsigned int *state, float mean) {
    if (mean <= 0) return 1;

    float p = 1.0f / (1.0f + mean);
    float u = rng_float(state);

    // Inverse transform: floor(log(u) / log(1-p))
    int result = (int)(logf(u + 1e-10f) / logf(1.0f - p + 1e-10f));
    return result < 1 ? 1 : result;
}

// Power law (Pareto) sampling
// P(X >= x) ~ x^(-alpha), heavy tailed
// Mean parameter controls the scale (x_min), alpha controls tail heaviness
// Small events common, but occasional large events occur
static int sample_power_law(unsigned int *state, float mean, float alpha) {
    if (mean <= 0) return 1;
    if (alpha <= 1.0f) alpha = 1.1f;  // Must be > 1 for finite mean

    // For Pareto: mean = alpha * x_min / (alpha - 1) when alpha > 1
    // So x_min = mean * (alpha - 1) / alpha
    float x_min = mean * (alpha - 1.0f) / alpha;
    if (x_min < 1.0f) x_min = 1.0f;

    float u = rng_float(state);

    // Inverse transform: x = x_min * (1 - u)^(-1/alpha)
    int result = (int)(x_min * powf(1.0f - u + 1e-10f, -1.0f / alpha));
    return result < 1 ? 1 : result;
}

// ============================================================================
// Distribution dispatch helpers
// ============================================================================

static int sample_count(unsigned int *state, CountDistribution dist, float mean, float dispersion) {
    switch (dist) {
        case DIST_NEGATIVE_BINOMIAL:
            return sample_negative_binomial(state, mean, dispersion);
        case DIST_POISSON:
        default:
            return sample_poisson(state, mean);
    }
}

static int sample_size(unsigned int *state, SizeDistribution dist, float mean, float power_law_alpha) {
    switch (dist) {
        case SIZE_GEOMETRIC:
            return sample_geometric(state, mean);
        case SIZE_POWER_LAW:
            return sample_power_law(state, mean, power_law_alpha);
        case SIZE_POISSON:
        default: {
            int s = sample_poisson(state, mean);
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

// Replace units[start:end) with `n_new` fresh units carved from `seq` (n_new *
// REPEAT_SIZE bytes). Frees the old units, shifts the tail, inserts the new ones.
// `seq` may be NULL when n_new == 0 (pure deletion).
static void replace_units_range(RepeatArray *arr, int start, int end,
                                const char *seq, int n_new) {
    int n_old = end - start;
    int delta = n_new - n_old;

    for (int i = start; i < end; i++) free(arr->units[i]);

    if (delta > 0) array_ensure_capacity(arr, arr->num_units + delta);

    // Shift the tail [end, num_units) to [start + n_new, ...)
    memmove(&arr->units[start + n_new], &arr->units[end],
            (arr->num_units - end) * sizeof(char *));

    for (int i = 0; i < n_new; i++) {
        arr->units[start + i] = alloc_unit(seq + (size_t)i * REPEAT_SIZE);
    }

    arr->num_units += delta;
}

// Copy base-pair range [a, b) out of the unit array into dst. Returns bytes written.
static long gather_bytes(char *dst, char **units, long a, long b) {
    long w = 0;
    while (a < b) {
        int u = (int)(a / REPEAT_SIZE);
        int off = (int)(a % REPEAT_SIZE);
        int n = REPEAT_SIZE - off;
        if ((long)n > b - a) n = (int)(b - a);
        memcpy(dst + w, units[u] + off, n);
        w += n; a += n;
    }
    return w;
}

// Tandem duplication of the base-pair window [char_start, char_end) (size a
// multiple of REPEAT_SIZE). char_start may be mid-unit, so re-splitting the
// rebuilt region on the REPEAT_SIZE grid produces chimeric units at the junctions
// -- mirrors censim's duplicate_at_position (the source of recombination diversity).
static void duplicate_at_position(RepeatArray *arr, long char_start, long char_end) {
    const int rs = REPEAT_SIZE;
    int unit_start = (int)(char_start / rs);
    int unit_end   = (int)((char_end + rs - 1) / rs);   // ceiling: last affected unit + 1
    long region_lo = (long)unit_start * rs;
    long region_hi = (long)unit_end * rs;
    long new_len = (region_hi - region_lo) + (char_end - char_start);  // multiple of rs

    char *buf = (char *)malloc(new_len);
    long w = 0;
    w += gather_bytes(buf + w, arr->units, region_lo, char_end);   // before + original segment
    w += gather_bytes(buf + w, arr->units, char_start, char_end);  // tandem copy
    w += gather_bytes(buf + w, arr->units, char_end, region_hi);   // after

    replace_units_range(arr, unit_start, unit_end, buf, (int)(new_len / rs));
    free(buf);
}

// Delete the base-pair window [char_start, char_end) (size a multiple of
// REPEAT_SIZE). Mid-unit endpoints merge into one chimeric unit. Mirrors censim's
// delete_at_position.
static void delete_at_position(RepeatArray *arr, long char_start, long char_end) {
    const int rs = REPEAT_SIZE;
    int unit_start = (int)(char_start / rs);
    int unit_end   = (int)((char_end + rs - 1) / rs);   // ceiling
    long region_lo = (long)unit_start * rs;
    long region_hi = (long)unit_end * rs;
    long new_len = (region_hi - region_lo) - (char_end - char_start);  // multiple of rs

    if (new_len == 0) {
        replace_units_range(arr, unit_start, unit_end, NULL, 0);
        return;
    }
    char *buf = (char *)malloc(new_len);
    long w = 0;
    w += gather_bytes(buf + w, arr->units, region_lo, char_start);  // kept head
    w += gather_bytes(buf + w, arr->units, char_end, region_hi);    // kept tail
    replace_units_range(arr, unit_start, unit_end, buf, (int)(new_len / rs));
    free(buf);
}

// ============================================================================
// Mutation functions
// ============================================================================

static void apply_snps(Simulation *sim) {
    unsigned int *rng = &sim->rng_state;
    int n_snps = sample_count(rng, sim->params.count_dist, sim->params.snp_rate, sim->params.nb_dispersion);

    for (int i = 0; i < n_snps; i++) {
        if (sim->array.num_units == 0) break;

        // Pick random unit and position
        int unit_idx = rng_int(rng, sim->array.num_units);
        int pos = rng_int(rng, REPEAT_SIZE);

        char *unit = sim->array.units[unit_idx];
        char old_base = unit[pos];

        // Pick a different base
        char new_base;
        do {
            new_base = BASES[rng_int(rng, 4)];
        } while (new_base == old_base);

        unit[pos] = new_base;
        sim->stats.snp_count++;
    }
}

static void apply_indels(Simulation *sim) {
    unsigned int *rng = &sim->rng_state;
    int n_indels = sample_count(rng, sim->params.count_dist, sim->params.indel_rate, sim->params.nb_dispersion);

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
        bool is_dup = rng_float(rng) < dup_prob;

        // Split the mean event size between dups and dels by the size ratio r, kept
        // log-symmetric around the central lambda: dup ~ lambda*sqrt(r), del ~
        // lambda/sqrt(r). r == 1 leaves both at lambda (size-symmetric, the default).
        float size_mean = sim->params.indel_size_lambda;
        float ratio = sim->params.dup_del_size_ratio;
        if (ratio > 0.0f && ratio != 1.0f) {
            float sr = sqrtf(ratio);
            size_mean = is_dup ? size_mean * sr : size_mean / sr;
        }

        // Sample size in repeat units (>=1), then convert to base pairs. The event
        // size is a whole multiple of REPEAT_SIZE so the frame is preserved, but the
        // START is an arbitrary base position (mid-unit), which produces chimeric
        // units at the junctions -- this is the dominant source of sequence
        // diversity, matching censim's arbitrary-position indels.
        int indel_units = sample_size(rng, sim->params.size_dist, size_mean, sim->params.power_law_alpha);

        long total_bp = (long)sim->array.num_units * REPEAT_SIZE;
        long char_start = rng_int(rng, (int)total_bp);
        long char_end = char_start + (long)indel_units * REPEAT_SIZE;

        // Out of bounds: can't place the event (array too small for this size). Skip.
        if (char_end >= total_bp) {
            continue;
        }

        if (is_dup) {
            // Check max size
            if (sim->params.bounding_enabled &&
                sim->array.num_units + indel_units > sim->params.max_array_size) {
                continue;
            }
            duplicate_at_position(&sim->array, char_start, char_end);
            sim->stats.dup_count++;
        } else {
            // Check min size
            if (sim->params.bounding_enabled &&
                sim->array.num_units - indel_units < sim->params.min_array_size) {
                continue;
            }
            delete_at_position(&sim->array, char_start, char_end);
            sim->stats.del_count++;
        }
    }

    // Check for collapse (independent of hard bounds; with hard bounds enabled the
    // array can't drop below min_array_size so this won't fire, matching the paper's
    // unbounded drift-to-collapse behaviour when bounding is off).
    if (sim->array.num_units < sim->params.collapse_threshold) {
        sim->stats.collapsed = true;
    }
}

// ============================================================================
// Public API
// ============================================================================

void sim_init(Simulation *sim, int initial_size, unsigned int seed) {
    // Seed per-trajectory RNG. Avoid 0 (xorshift fixed point).
    sim->rng_state = seed ? seed : 0x9E3779B9u;

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
    sim->params.dup_del_size_ratio = DEFAULT_DUP_DEL_SIZE_RATIO;
    sim->params.snp_rate = DEFAULT_SNP_RATE;
    sim->params.min_array_size = DEFAULT_MIN_ARRAY_SIZE;
    sim->params.max_array_size = DEFAULT_MAX_ARRAY_SIZE;
    sim->params.bounding_enabled = true;
    sim->params.collapse_threshold = DEFAULT_COLLAPSE_THRESHOLD;
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

#include "colorizer.h"
#include "config.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ENCODING_SIZE (REPEAT_SIZE * 4)  // 712
#define CACHE_INITIAL_SIZE 4096

// ============================================================================
// Random number generation (same as simulation but separate state)
// ============================================================================

static unsigned int col_rng_state = 42;

static void col_rng_seed(unsigned int seed) {
    col_rng_state = seed;
}

static float col_rng_gauss(void) {
    // Box-Muller transform for Gaussian random numbers
    static int has_spare = 0;
    static float spare;

    if (has_spare) {
        has_spare = 0;
        return spare;
    }

    float u, v, s;
    do {
        col_rng_state ^= col_rng_state << 13;
        col_rng_state ^= col_rng_state >> 17;
        col_rng_state ^= col_rng_state << 5;
        u = 2.0f * ((float)(col_rng_state & 0x7FFFFFFF) / (float)0x7FFFFFFF) - 1.0f;

        col_rng_state ^= col_rng_state << 13;
        col_rng_state ^= col_rng_state >> 17;
        col_rng_state ^= col_rng_state << 5;
        v = 2.0f * ((float)(col_rng_state & 0x7FFFFFFF) / (float)0x7FFFFFFF) - 1.0f;

        s = u * u + v * v;
    } while (s >= 1.0f || s == 0.0f);

    float mul = sqrtf(-2.0f * logf(s) / s);
    spare = v * mul;
    has_spare = 1;
    return u * mul;
}

// ============================================================================
// Orthogonalization (Gram-Schmidt)
// ============================================================================

static void orthogonalize(float matrix[][3], int rows) {
    for (int col = 0; col < 3; col++) {
        // Subtract projections onto previous columns
        for (int prev = 0; prev < col; prev++) {
            float dot = 0.0f, norm_sq = 0.0f;
            for (int row = 0; row < rows; row++) {
                dot += matrix[row][col] * matrix[row][prev];
                norm_sq += matrix[row][prev] * matrix[row][prev];
            }
            if (norm_sq > 0) {
                float scale = dot / norm_sq;
                for (int row = 0; row < rows; row++) {
                    matrix[row][col] -= scale * matrix[row][prev];
                }
            }
        }

        // Normalize column
        float norm = 0.0f;
        for (int row = 0; row < rows; row++) {
            norm += matrix[row][col] * matrix[row][col];
        }
        norm = sqrtf(norm);
        if (norm > 0) {
            for (int row = 0; row < rows; row++) {
                matrix[row][col] /= norm;
            }
        }
    }
}

// ============================================================================
// One-hot encoding
// ============================================================================

static void one_hot_encode(const char *seq, float *out) {
    memset(out, 0, ENCODING_SIZE * sizeof(float));

    for (int i = 0; i < REPEAT_SIZE; i++) {
        int base_idx;
        switch (seq[i]) {
            case 'A': base_idx = 0; break;
            case 'C': base_idx = 1; break;
            case 'G': base_idx = 2; break;
            case 'T': base_idx = 3; break;
            default: continue;
        }
        out[i * 4 + base_idx] = 1.0f;
    }
}

// ============================================================================
// Hashing for cache
// ============================================================================

static unsigned int hash_sequence(const char *seq) {
    unsigned int hash = 2166136261u;
    for (int i = 0; i < REPEAT_SIZE; i++) {
        hash ^= (unsigned char)seq[i];
        hash *= 16777619u;
    }
    return hash;
}

// ============================================================================
// Fixed bounds computation
// ============================================================================

static void compute_fixed_bounds(Colorizer *c) {
    // Sample random sequences to estimate the range
    int n_samples = 1000;
    float *raw = (float *)malloc(n_samples * 3 * sizeof(float));
    float encoding[ENCODING_SIZE];

    col_rng_seed(col_rng_state + 1000);  // Different seed for sampling

    for (int s = 0; s < n_samples; s++) {
        // Generate random sequence
        char seq[REPEAT_SIZE + 1];
        for (int i = 0; i < REPEAT_SIZE; i++) {
            col_rng_state ^= col_rng_state << 13;
            col_rng_state ^= col_rng_state >> 17;
            col_rng_state ^= col_rng_state << 5;
            seq[i] = BASES[col_rng_state % 4];
        }
        seq[REPEAT_SIZE] = '\0';

        // Encode and project
        one_hot_encode(seq, encoding);
        for (int ch = 0; ch < 3; ch++) {
            float val = 0.0f;
            for (int j = 0; j < ENCODING_SIZE; j++) {
                val += encoding[j] * c->projection[j][ch];
            }
            raw[s * 3 + ch] = val;
        }
    }

    // Sort each channel and use percentiles
    for (int ch = 0; ch < 3; ch++) {
        float channel_vals[1000];
        for (int s = 0; s < n_samples; s++) {
            channel_vals[s] = raw[s * 3 + ch];
        }

        // Simple bubble sort (only 1000 elements)
        for (int i = 0; i < n_samples - 1; i++) {
            for (int j = 0; j < n_samples - i - 1; j++) {
                if (channel_vals[j] > channel_vals[j + 1]) {
                    float tmp = channel_vals[j];
                    channel_vals[j] = channel_vals[j + 1];
                    channel_vals[j + 1] = tmp;
                }
            }
        }

        // 1st and 99th percentile with padding
        c->min_vals[ch] = channel_vals[10] - 0.5f;
        c->max_vals[ch] = channel_vals[989] + 0.5f;
    }

    free(raw);
}

// ============================================================================
// Public API
// ============================================================================

void colorizer_init(Colorizer *c, unsigned int seed) {
    col_rng_seed(seed);

    // Generate random Gaussian matrix
    for (int row = 0; row < ENCODING_SIZE; row++) {
        for (int col = 0; col < 3; col++) {
            c->projection[row][col] = col_rng_gauss();
        }
    }

    // Orthogonalize
    orthogonalize(c->projection, ENCODING_SIZE);

    // Compute fixed bounds
    compute_fixed_bounds(c);

    // Initialize cache
    c->cache_capacity = CACHE_INITIAL_SIZE;
    c->cache_size = 0;
    c->cache_hashes = (unsigned int *)calloc(c->cache_capacity, sizeof(unsigned int));
    c->cache_colors = (Color *)malloc(c->cache_capacity * sizeof(Color));
}

void colorizer_free(Colorizer *c) {
    free(c->cache_hashes);
    free(c->cache_colors);
    c->cache_hashes = NULL;
    c->cache_colors = NULL;
    c->cache_size = 0;
    c->cache_capacity = 0;
}

void colorizer_clear_cache(Colorizer *c) {
    memset(c->cache_hashes, 0, c->cache_capacity * sizeof(unsigned int));
    c->cache_size = 0;
}

Color colorizer_get_color(Colorizer *c, const char *seq) {
    unsigned int hash = hash_sequence(seq);

    // Check cache
    unsigned int idx = hash % c->cache_capacity;
    int probes = 0;
    while (c->cache_hashes[idx] != 0 && probes < c->cache_capacity) {
        if (c->cache_hashes[idx] == hash) {
            return c->cache_colors[idx];
        }
        idx = (idx + 1) % c->cache_capacity;
        probes++;
    }

    // Not in cache - compute color
    float encoding[ENCODING_SIZE];
    one_hot_encode(seq, encoding);

    float rgb[3];
    for (int ch = 0; ch < 3; ch++) {
        float val = 0.0f;
        for (int j = 0; j < ENCODING_SIZE; j++) {
            val += encoding[j] * c->projection[j][ch];
        }

        // Normalize to [0, 1]
        float range = c->max_vals[ch] - c->min_vals[ch];
        if (range > 0) {
            val = (val - c->min_vals[ch]) / range;
        }
        if (val < 0) val = 0;
        if (val > 1) val = 1;
        rgb[ch] = val;
    }

    Color color = {
        (unsigned char)(rgb[0] * 255),
        (unsigned char)(rgb[1] * 255),
        (unsigned char)(rgb[2] * 255),
        255
    };

    // Add to cache if there's room
    if (c->cache_size < c->cache_capacity * 3 / 4) {
        idx = hash % c->cache_capacity;
        while (c->cache_hashes[idx] != 0) {
            idx = (idx + 1) % c->cache_capacity;
        }
        c->cache_hashes[idx] = hash;
        c->cache_colors[idx] = color;
        c->cache_size++;
    }

    return color;
}

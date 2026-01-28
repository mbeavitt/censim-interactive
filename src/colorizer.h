#ifndef COLORIZER_H
#define COLORIZER_H

#include <raylib.h>

// Orthogonal projection colorizer
typedef struct {
    float projection[712][3];  // 178*4 x 3 orthogonal matrix
    float min_vals[3];         // Fixed normalization bounds
    float max_vals[3];

    // Color cache (simple hash table)
    unsigned int *cache_hashes;
    Color *cache_colors;
    int cache_size;
    int cache_capacity;
} Colorizer;

void colorizer_init(Colorizer *c, unsigned int seed);
void colorizer_free(Colorizer *c);
Color colorizer_get_color(Colorizer *c, const char *seq);
void colorizer_clear_cache(Colorizer *c);

#endif // COLORIZER_H

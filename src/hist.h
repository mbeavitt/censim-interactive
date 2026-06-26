#ifndef HIST_H
#define HIST_H

// Minimal fixed-range histogram, reused for both per-array metric distributions
// (one value per trajectory) and per-HOR metric distributions (many values per
// trajectory folded into shared bins). Supports linear or log10 binning.

typedef struct {
    float min, max;       // value range (for log scale these are raw values, > 0)
    int   nbins;
    int   log_scale;      // 1 = bin in log10 space, 0 = linear
    long *counts;         // nbins entries
    long  total;          // values added (including under/overflow)
    long  underflow;      // values < min
    long  overflow;       // values > max
} Histogram;

void  hist_init(Histogram *h, float min, float max, int nbins, int log_scale);
void  hist_free(Histogram *h);
void  hist_reset(Histogram *h);
void  hist_add(Histogram *h, float value);
void  hist_add_weighted(Histogram *h, float value, long weight);
void  hist_merge(Histogram *dst, const Histogram *src);  // sum counts; bins must match
float hist_bin_center(const Histogram *h, int i);        // value at center of bin i
float hist_bin_lo(const Histogram *h, int i);            // left edge of bin i
long  hist_max_count(const Histogram *h);                // tallest bin (for plot scaling)

#endif // HIST_H

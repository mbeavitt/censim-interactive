#include "hist.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void hist_init(Histogram *h, float min, float max, int nbins, int log_scale) {
    h->min = min;
    h->max = max;
    h->nbins = nbins;
    h->log_scale = log_scale;
    h->counts = (long *)calloc(nbins, sizeof(long));
    h->total = h->underflow = h->overflow = 0;
}

void hist_free(Histogram *h) {
    free(h->counts);
    h->counts = NULL;
    h->nbins = 0;
}

void hist_reset(Histogram *h) {
    if (h->counts) memset(h->counts, 0, h->nbins * sizeof(long));
    h->total = h->underflow = h->overflow = 0;
}

// Map a value to a bin index, or -1 (underflow) / nbins (overflow).
static int hist_bin_index(const Histogram *h, float value) {
    float lo = h->min, hi = h->max, v = value;
    if (h->log_scale) {
        // A log axis can't place values <= 0 (e.g. block gap == 0 for overlapping
        // HORs, which are very common). Fold those -- and anything below min -- into
        // the first bin so they show as a left-edge spike instead of vanishing into
        // underflow. Only genuine overflow (> max) is flagged.
        if (v <= 0.0f) return 0;
        lo = log10f(h->min);
        hi = log10f(h->max);
        v  = log10f(v);
        if (v < lo) return 0;
        if (v > hi) return h->nbins;
        float frac = (v - lo) / (hi - lo);
        int idx = (int)(frac * h->nbins);
        if (idx >= h->nbins) idx = h->nbins - 1;
        if (idx < 0) idx = 0;
        return idx;
    }
    if (v < lo) return -1;
    if (v > hi) return h->nbins;
    float frac = (v - lo) / (hi - lo);
    int idx = (int)(frac * h->nbins);
    if (idx >= h->nbins) idx = h->nbins - 1;  // include the upper edge
    if (idx < 0) idx = 0;
    return idx;
}

void hist_add_weighted(Histogram *h, float value, long weight) {
    int idx = hist_bin_index(h, value);
    h->total += weight;
    if (idx < 0)            { h->underflow += weight; return; }
    if (idx >= h->nbins)    { h->overflow  += weight; return; }
    h->counts[idx] += weight;
}

void hist_add(Histogram *h, float value) {
    hist_add_weighted(h, value, 1);
}

void hist_merge(Histogram *dst, const Histogram *src) {
    if (dst->nbins != src->nbins) return;
    for (int i = 0; i < dst->nbins; i++) dst->counts[i] += src->counts[i];
    dst->total     += src->total;
    dst->underflow += src->underflow;
    dst->overflow  += src->overflow;
}

float hist_bin_lo(const Histogram *h, int i) {
    if (h->log_scale) {
        float lo = log10f(h->min), hi = log10f(h->max);
        return powf(10.0f, lo + (hi - lo) * i / h->nbins);
    }
    return h->min + (h->max - h->min) * i / h->nbins;
}

float hist_bin_center(const Histogram *h, int i) {
    if (h->log_scale) {
        // geometric center of the bin
        return sqrtf(hist_bin_lo(h, i) * hist_bin_lo(h, i + 1));
    }
    return 0.5f * (hist_bin_lo(h, i) + hist_bin_lo(h, i + 1));
}

long hist_max_count(const Histogram *h) {
    long m = 0;
    for (int i = 0; i < h->nbins; i++) if (h->counts[i] > m) m = h->counts[i];
    return m;
}

#ifndef HOR_H
#define HOR_H

#include "simulation.h"
#include "hist.h"

// In-engine higher-order-repeat (HOR) detection.
//
// This is a port of Piotr's TRASH HORT logic, specialised for simulated arrays.
// Because every simulated unit is exactly REPEAT_SIZE bp (SNPs substitute in
// place; dup/del move whole units), there are no length-varying repeats and
// therefore no need for a multiple-sequence alignment / mafft step: comparing
// two units is a plain capped Hamming distance over REPEAT_SIZE columns.
//
// The sim never inverts, so only direct (parallel) repeats exist; detection is a
// single upper-triangular diagonal scan of the self-similarity matrix. A HOR is a
// maximal run of consecutive units (i, i+1, ...) that stays pairwise-similar to
// another run (i+d, i+1+d, ...), of length >= HOR_CUTOFF, with every aligned pair
// within HOR_THRESHOLD substitutions.
//
// Detection parameters match the report / TRASH flags (--hor_threshold=3,
// --min_hor_len=3).

// Detection parameters. Overridable at compile time (-D) so we can match the
// ground-truth HORT run, whose threshold is percentage-scaled:
// floor(hor_threshold% * median_width / 100) = floor(4 * 178 / 100) = 7.
#ifndef HOR_THRESHOLD
#define HOR_THRESHOLD 7   // max substitutions allowed between an aligned pair
#endif
#ifndef HOR_CUTOFF
#define HOR_CUTOFF    3   // minimum block length (units) to count as a HOR
#endif

// One detected HOR (a pair of similar blocks A and B).
typedef struct {
    int    block_size;    // units per block (the diagonal run length)
    int    block_gap;     // units between end of A and start of B (0 if overlapping)
    float  similarity;    // 1 / (1 + mean per-position substitutions between A and B)
    float  diversity;     // unique units in block A / block_size  (1/n .. 1)
    double composite;     // block_gap * similarity * block_size * diversity
} HorBlock;

// Per-array summary (one record per trajectory).
typedef struct {
    long   num_hors;
    int    num_units;
    double array_kb;       // num_units * REPEAT_SIZE / 1000
    double hors_per_kb;
    double mean_size;
    double mean_gap;
    double mean_similarity;
    double mean_diversity;
    double mean_composite;
    long   diversity_eq_one;   // # HORs whose block A is all-unique (diversity == 1)
} HorStats;

// Optional per-HOR distributions, folded in during the scan. Any field may be
// NULL to skip it. Callers create these with the bin ranges they want to plot.
typedef struct {
    Histogram *size;        // block size (units), typically log-scaled
    Histogram *gap;         // block gap (units), typically log-scaled
    Histogram *similarity;  // 0..1 linear
    Histogram *diversity;   // 0..1 linear
    Histogram *composite;   // log-scaled
} HorHistSet;

// Scan an array, filling `stats` (required) and folding per-HOR values into
// `hists` (may be NULL).
void hor_scan(const RepeatArray *array, HorStats *stats, HorHistSet *hists);

// Scan and return every detected HOR in a freshly malloc'd array (caller frees).
// Intended for validation / inspection on small arrays. *out_count is set.
HorBlock *hor_scan_collect(const RepeatArray *array, long *out_count);

#endif // HOR_H

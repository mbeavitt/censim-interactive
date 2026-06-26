#include "hor.h"
#include "config.h"
#include <stdlib.h>
#include <string.h>
// Substitution count between two aligned units (Piotr's compareAB). Exact, plain
// byte loop; optimize later if needed.
static int hamming(const char *a, const char *b) {
    int d = 0;
    for (int i = 0; i < REPEAT_SIZE; i++) {
        d += (a[i] != b[i]);
    }
    return d;
}

// FNV-1a over a unit; used to count unique sequences within a block cheaply.
static unsigned int unit_hash(const char *seq) {
    unsigned int h = 2166136261u;
    for (int i = 0; i < REPEAT_SIZE; i++) {
        h ^= (unsigned char)seq[i];
        h *= 16777619u;
    }
    return h;
}

static int cmp_uint(const void *a, const void *b) {
    unsigned int x = *(const unsigned int *)a, y = *(const unsigned int *)b;
    return (x > y) - (x < y);
}

// Count distinct units in block [start, start+len). Hash + sort + dedup.
// Hash collisions could in principle merge two distinct units, but with a 32-bit
// FNV hash over small blocks this is negligible for our purposes.
static int count_unique_block(char **units, int start, int len) {
    if (len <= 1) return len;
    unsigned int stackbuf[256];
    unsigned int *h = (len <= 256) ? stackbuf
                                   : (unsigned int *)malloc(len * sizeof(unsigned int));
    for (int i = 0; i < len; i++) h[i] = unit_hash(units[start + i]);
    qsort(h, len, sizeof(unsigned int), cmp_uint);
    int unique = 1;
    for (int i = 1; i < len; i++) if (h[i] != h[i - 1]) unique++;
    if (h != stackbuf) free(h);
    return unique;
}

// Build a HorBlock from a detected run: block A starts at aStart, length `run`,
// block B is the same run shifted by diagonal `d`. `run_snv` is the summed
// substitution count over the run's aligned pairs.
static HorBlock make_block(char **units, int aStart, int run, int d, long run_snv) {
    HorBlock b;
    b.block_size = run;
    int gap = d - run;                 // start_B - end_A - 1
    b.block_gap = gap < 0 ? 0 : gap;   // overlapping blocks -> gap 0 (per report)
    float mean_snv = (float)run_snv / (float)run;
    b.similarity = 1.0f / (1.0f + mean_snv);
    int uniq = count_unique_block(units, aStart, run);
    b.diversity = (float)uniq / (float)run;
    b.composite = (double)b.block_gap * b.similarity * b.block_size * b.diversity;
    return b;
}

static void stats_accumulate(HorStats *s, const HorBlock *b) {
    s->num_hors++;
    s->mean_size       += b->block_size;
    s->mean_gap        += b->block_gap;
    s->mean_similarity += b->similarity;
    s->mean_diversity  += b->diversity;
    s->mean_composite  += b->composite;
    if (b->diversity >= 1.0f) s->diversity_eq_one++;
}

static void hists_fold(HorHistSet *h, const HorBlock *b) {
    if (!h) return;
    if (h->size)       hist_add(h->size, (float)b->block_size);
    if (h->gap)        hist_add(h->gap, (float)b->block_gap);
    if (h->similarity) hist_add(h->similarity, b->similarity);
    if (h->diversity)  hist_add(h->diversity, b->diversity);
    if (h->composite)  hist_add(h->composite, (float)b->composite);
}

void hor_scan(const RepeatArray *array, HorStats *stats, HorHistSet *hists) {
    memset(stats, 0, sizeof(*stats));
    int n = array->num_units;
    stats->num_units = n;
    stats->array_kb = (double)n * REPEAT_SIZE / 1000.0;
    if (n <= HOR_CUTOFF) return;

    char **u = array->units;

    // Walk each diagonal offset d; find maximal pairwise-similar runs along it.
    for (int d = 1; d < n; d++) {
        int  run = 0;
        long run_snv = 0;
        int  limit = n - d;
        for (int i = 0; i < limit; i++) {
            int hd = hamming(u[i], u[i + d]);
            if (hd <= HOR_THRESHOLD) {
                run++;
                run_snv += hd;
            } else {
                if (run >= HOR_CUTOFF) {
                    HorBlock b = make_block(u, i - run, run, d, run_snv);
                    stats_accumulate(stats, &b);
                    hists_fold(hists, &b);
                }
                run = 0;
                run_snv = 0;
            }
        }
        if (run >= HOR_CUTOFF) {  // run reached the diagonal's edge still open
            HorBlock b = make_block(u, limit - run, run, d, run_snv);
            stats_accumulate(stats, &b);
            hists_fold(hists, &b);
        }
    }

    if (stats->num_hors > 0) {
        double inv = 1.0 / (double)stats->num_hors;
        stats->mean_size       *= inv;
        stats->mean_gap        *= inv;
        stats->mean_similarity *= inv;
        stats->mean_diversity  *= inv;
        stats->mean_composite  *= inv;
    }
    stats->hors_per_kb = (stats->array_kb > 0.0)
                       ? (double)stats->num_hors / stats->array_kb : 0.0;
}

HorBlock *hor_scan_collect(const RepeatArray *array, long *out_count) {
    int n = array->num_units;
    long cap = 1024, count = 0;
    HorBlock *out = (HorBlock *)malloc(cap * sizeof(HorBlock));
    char **u = array->units;

    for (int d = 1; d < n; d++) {
        int  run = 0;
        long run_snv = 0;
        int  limit = n - d;
        for (int i = 0; i < limit; i++) {
            int hd = hamming(u[i], u[i + d]);
            if (hd <= HOR_THRESHOLD) {
                run++;
                run_snv += hd;
            } else {
                if (run >= HOR_CUTOFF) {
                    if (count == cap) { cap *= 2; out = realloc(out, cap * sizeof(HorBlock)); }
                    out[count++] = make_block(u, i - run, run, d, run_snv);
                }
                run = 0;
                run_snv = 0;
            }
        }
        if (run >= HOR_CUTOFF) {
            if (count == cap) { cap *= 2; out = realloc(out, cap * sizeof(HorBlock)); }
            out[count++] = make_block(u, limit - run, run, d, run_snv);
        }
    }

    *out_count = count;
    return out;
}

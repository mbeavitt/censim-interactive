#include "hor.h"
#include "config.h"
#include <stdlib.h>
#include <string.h>
// Substitution count between two aligned units (Piotr's compareAB).
// Portable scalar version (also the fallback on non-x86 / pre-AVX2 CPUs).
static int hamming_scalar(const char *a, const char *b) {
    int d = 0;
    for (int i = 0; i < REPEAT_SIZE; i++) {
        d += (a[i] != b[i]);
    }
    return d;
}

#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
// AVX2 base-mismatch count: compare 32 bytes at a time (cmpeq_epi8 -> per-byte
// equality mask), count the unequal bytes via popcount of the inverted mask.
// Compiled with AVX2 codegen via the target attribute regardless of global flags;
// only invoked when the CPU supports it (see hamming()), so the binary stays
// portable. Same idea as the hamming256 popcount trick in kmer-variance.
__attribute__((target("avx2")))
static int hamming_avx2(const char *a, const char *b) {
    int diff = 0, i = 0;
    for (; i + 32 <= REPEAT_SIZE; i += 32) {
        __m256i va = _mm256_loadu_si256((const __m256i *)(a + i));
        __m256i vb = _mm256_loadu_si256((const __m256i *)(b + i));
        unsigned int eq = (unsigned int)_mm256_movemask_epi8(_mm256_cmpeq_epi8(va, vb));
        diff += 32 - __builtin_popcount(eq);   // bits set in eq = equal bytes
    }
    for (; i < REPEAT_SIZE; i++) diff += (a[i] != b[i]);
    return diff;
}

static int hamming(const char *a, const char *b) {
    static int use_avx2 = -1;
    if (use_avx2 < 0) use_avx2 = __builtin_cpu_supports("avx2");
    return use_avx2 ? hamming_avx2(a, b) : hamming_scalar(a, b);
}
#else
static int hamming(const char *a, const char *b) { return hamming_scalar(a, b); }
#endif

// FNV-1a over a unit; used to count unique sequences within a block cheaply.
static unsigned int unit_hash(const char *seq) {
    unsigned int h = 2166136261u;
    for (int i = 0; i < REPEAT_SIZE; i++) {
        h ^= (unsigned char)seq[i];
        h *= 16777619u;
    }
    return h;
}

// Count distinct units in block [start, start+len) using precomputed per-unit
// hashes, via an open-addressing hash set with linear probing: O(len), no sort.
// 0 is the empty slot; an actual hash of 0 is tracked separately. Hash collisions
// could in principle merge two distinct units, but with a 32-bit FNV hash over
// small blocks this is negligible.
static int count_unique_block(const unsigned int *hashes, int start, int len) {
    if (len <= 1) return len;
    int cap = 16;
    while (cap < len * 2) cap <<= 1;          // power-of-two, load factor < 0.5
    unsigned int stackbuf[2048];              // covers len up to 1024 without malloc
    unsigned int *tab = (cap <= 2048) ? stackbuf
                                      : (unsigned int *)malloc((size_t)cap * sizeof(unsigned int));
    memset(tab, 0, (size_t)cap * sizeof(unsigned int));
    unsigned int mask = (unsigned int)cap - 1u;
    int unique = 0, have_zero = 0;
    for (int i = 0; i < len; i++) {
        unsigned int hv = hashes[start + i];
        if (hv == 0) { if (!have_zero) { have_zero = 1; unique++; } continue; }
        unsigned int idx = hv & mask;
        while (tab[idx] != 0 && tab[idx] != hv) idx = (idx + 1) & mask;
        if (tab[idx] == 0) { tab[idx] = hv; unique++; }
    }
    if (tab != stackbuf) free(tab);
    return unique;
}

// Build a HorBlock from a detected run: block A starts at aStart, length `run`,
// block B is the same run shifted by diagonal `d`. `run_snv` is the summed
// substitution count over the run's aligned pairs. `hashes` is the array-wide
// table of per-unit hashes (precomputed once per scan).
static HorBlock make_block(const unsigned int *hashes, int aStart, int run, int d, long run_snv) {
    HorBlock b;
    b.block_size = run;
    int gap = d - run;                 // start_B - end_A - 1
    b.block_gap = gap < 0 ? 0 : gap;   // overlapping blocks -> gap 0 (per report)
    float mean_snv = (float)run_snv / (float)run;
    b.similarity = 1.0f / (1.0f + mean_snv);
    int uniq = count_unique_block(hashes, aStart, run);
    b.diversity = (float)uniq / (float)run;
    b.composite = (double)b.block_gap * b.similarity * b.block_size * b.diversity;
    return b;
}

// Precompute the FNV hash of every unit once (the per-HOR diversity calc reuses
// these instead of re-hashing the same units hundreds of millions of times).
static unsigned int *precompute_hashes(char **units, int n) {
    unsigned int *hashes = (unsigned int *)malloc((size_t)n * sizeof(unsigned int));
    for (int i = 0; i < n; i++) hashes[i] = unit_hash(units[i]);
    return hashes;
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

void hor_scan(const RepeatArray *array, HorStats *stats, HorHistSet *hists,
              const volatile bool *cancel) {
    memset(stats, 0, sizeof(*stats));
    int n = array->num_units;
    stats->num_units = n;
    stats->array_kb = (double)n * REPEAT_SIZE / 1000.0;
    if (n <= HOR_CUTOFF) return;

    char **u = array->units;
    unsigned int *hashes = precompute_hashes(u, n);

    // Walk each diagonal offset d; find maximal pairwise-similar runs along it.
    for (int d = 1; d < n; d++) {
        if (cancel && *cancel) break;  // abort early on cancellation
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
                    HorBlock b = make_block(hashes, i - run, run, d, run_snv);
                    stats_accumulate(stats, &b);
                    hists_fold(hists, &b);
                }
                run = 0;
                run_snv = 0;
            }
        }
        if (run >= HOR_CUTOFF) {  // run reached the diagonal's edge still open
            HorBlock b = make_block(hashes, limit - run, run, d, run_snv);
            stats_accumulate(stats, &b);
            hists_fold(hists, &b);
        }
    }
    free(hashes);

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
    unsigned int *hashes = precompute_hashes(u, n);

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
                    out[count++] = make_block(hashes, i - run, run, d, run_snv);
                }
                run = 0;
                run_snv = 0;
            }
        }
        if (run >= HOR_CUTOFF) {
            if (count == cap) { cap *= 2; out = realloc(out, cap * sizeof(HorBlock)); }
            out[count++] = make_block(hashes, limit - run, run, d, run_snv);
        }
    }

    free(hashes);
    *out_count = count;
    return out;
}

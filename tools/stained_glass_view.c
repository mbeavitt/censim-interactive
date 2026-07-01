// stained_glass_view.c
// ---------------------------------------------------------------------------
// Standalone raylib viewer for "stained glass" self-identity plots of a
// censim FASTA export. Each repeat unit is one row/column; cell (i,j) is the
// estimated sequence identity between units i and j.
//
// The all-vs-all is O(N^2), so instead of aligning every pair we approximate
// identity from k-mer content using MinHash sketches:
//   - each unit -> s independent min-hashes over its 16-mers (packed 2 bits/base)
//   - Jaccard(i,j) ~ fraction of the s sketch slots that match
//   - identity via the Mash estimator: id = 1 + (1/k) * ln(2J/(1+J))
// Comparing two units is then s integer compares instead of an alignment.
//
// Build (via premake): premake5 gmake2 && make config=release stained_glass
// Or directly:
//   cc -O2 -std=gnu99 tools/stained_glass_view.c -o bin/stained_glass \
//      -lraylib -lGL -lm -lpthread -ldl -lrt -lX11
//
// Run: ./bin/stained_glass example.fasta
// Controls: mouse wheel = zoom, left-drag = pan, R = reset view, ESC = quit.
// ---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>

#include <raylib.h>

static int num_threads(void) {
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    if (n < 1) n = 1;
    if (n > 64) n = 64;
    return (int)n;
}

static int g_quiet = 0;   // silence per-unit progress prints in benchmark mode

static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1.0e6;
}

#define KMER_K      16          // 16 bases pack exactly into a uint32 (2 bits each)
#define SKETCH_S    64          // MinHash sketch size (identity-estimate resolution)
#define MAX_DIM     1400        // cap on the rendered matrix side; larger N is binned

// ----------------------------------------------------------------------------
// FASTA loading
// ----------------------------------------------------------------------------

typedef struct {
    char **seqs;
    int    count;
} Fasta;

static Fasta load_fasta(const char *path) {
    Fasta fa = {0};
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "Cannot open FASTA: %s\n", path); exit(1); }

    int cap = 1024;
    fa.seqs = malloc(sizeof(char *) * cap);

    char *cur = NULL;
    size_t cur_len = 0, cur_cap = 0;
    char line[8192];

    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '>') {
            if (cur) { fa.seqs[fa.count++] = cur; cur = NULL; cur_len = cur_cap = 0; }
            if (fa.count >= cap) { cap *= 2; fa.seqs = realloc(fa.seqs, sizeof(char *) * cap); }
        } else {
            size_t n = strcspn(line, "\r\n");
            if (cur_len + n + 1 > cur_cap) {
                cur_cap = (cur_len + n + 1) * 2;
                cur = realloc(cur, cur_cap);
            }
            memcpy(cur + cur_len, line, n);
            cur_len += n;
            cur[cur_len] = '\0';
        }
    }
    if (cur) fa.seqs[fa.count++] = cur;
    fclose(f);
    return fa;
}

// ----------------------------------------------------------------------------
// MinHash signatures
// ----------------------------------------------------------------------------

static inline uint32_t fmix32(uint32_t h) {
    h ^= h >> 16; h *= 0x85ebca6b;
    h ^= h >> 13; h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

static inline int base2bit(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return -1;   // N / ambiguous -> break the k-mer window
    }
}

// Sketch a single unit into s[0..SKETCH_S).
static void sketch_unit(const char *seq, const uint32_t *seeds, uint32_t *s) {
    for (int h = 0; h < SKETCH_S; h++) s[h] = 0xFFFFFFFFu;
    uint32_t code = 0;
    int valid = 0;                        // consecutive valid bases in the window
    for (const char *p = seq; *p; p++) {
        int b = base2bit(*p);
        if (b < 0) { valid = 0; code = 0; continue; }
        code = (code << 2) | (uint32_t)b;  // k=16 => full 32-bit window, no mask needed
        if (++valid >= KMER_K) {
            for (int h = 0; h < SKETCH_S; h++) {
                uint32_t hv = fmix32(code ^ seeds[h]);
                if (hv < s[h]) s[h] = hv;
            }
        }
    }
}

typedef struct {
    const Fasta   *fa;
    uint32_t      *sig;
    const uint32_t *seeds;
    int            tid, nthreads;
} SigArg;

static void *sig_worker(void *p) {
    SigArg *a = p;
    for (int u = a->tid; u < a->fa->count; u += a->nthreads)
        sketch_unit(a->fa->seqs[u], a->seeds, a->sig + (size_t)u * SKETCH_S);
    return NULL;
}

// One flat array of N*SKETCH_S uint32 min-hashes. Units are independent, so we
// simply fan them out across threads (disjoint output rows, no locking).
static uint32_t *build_signatures(const Fasta *fa) {
    uint32_t seeds[SKETCH_S];
    for (int h = 0; h < SKETCH_S; h++) seeds[h] = fmix32((uint32_t)h * 0x9E3779B1u + 1u);

    uint32_t *sig = malloc((size_t)fa->count * SKETCH_S * sizeof(uint32_t));

    int T = num_threads();
    if (T > fa->count) T = fa->count;
    if (T < 1) T = 1;

    pthread_t *th = malloc(sizeof(pthread_t) * T);
    SigArg    *args = malloc(sizeof(SigArg) * T);
    for (int t = 0; t < T; t++) {
        args[t] = (SigArg){ fa, sig, seeds, t, T };
        pthread_create(&th[t], NULL, sig_worker, &args[t]);
    }
    for (int t = 0; t < T; t++) pthread_join(th[t], NULL);
    free(th); free(args);
    return sig;
}

// ----------------------------------------------------------------------------
// Build the (possibly binned) identity matrix
// ----------------------------------------------------------------------------

// Cap on matrix threads: each keeps a private dim*dim accumulation buffer
// (~24 MB at dim=1400), so bound the total memory footprint.
#define MATRIX_MAX_THREADS 16

typedef struct {
    const uint32_t *sig;
    const int      *bin;
    const float    *lut;      // identity indexed by number of matching sketch slots
    int             N, dim;
    int             tid, nthreads;
    double         *sum;      // private dim*dim accumulator
    uint32_t       *cnt;      // private dim*dim counts
} MatArg;

static void *mat_worker(void *p) {
    MatArg *a = p;
    int N = a->N, dim = a->dim;
    // Interleave rows across threads so each gets a balanced mix of long (small
    // i) and short (large i) inner spans. Every unordered pair {i,j} is handled
    // exactly once, by the thread that owns row i.
    for (int i = a->tid; i < N; i += a->nthreads) {
        const uint32_t *si = a->sig + (size_t)i * SKETCH_S;
        int bi = a->bin[i];
        for (int j = i; j < N; j++) {
            const uint32_t *sj = a->sig + (size_t)j * SKETCH_S;
            int matches = 0;
            for (int h = 0; h < SKETCH_S; h++) matches += (si[h] == sj[h]);
            float id = a->lut[matches];
            int bj = a->bin[j];
            size_t x = (size_t)bi * dim + bj;
            size_t y = (size_t)bj * dim + bi;
            a->sum[x] += id; a->cnt[x]++;
            if (x != y) { a->sum[y] += id; a->cnt[y]++; }
        }
    }
    return NULL;
}

// Returns a dim*dim float matrix of average identities; sets *out_dim.
static float *build_matrix(const Fasta *fa, const uint32_t *sig, int *out_dim) {
    int N = fa->count;
    int dim = N < MAX_DIM ? N : MAX_DIM;
    *out_dim = dim;

    // Identity is a function of the integer match count only (0..SKETCH_S), so
    // precompute it once instead of calling logf per pair.
    float lut[SKETCH_S + 1];
    for (int m = 0; m <= SKETCH_S; m++) {
        if (m == 0)            { lut[m] = 0.0f; continue; }
        if (m == SKETCH_S)     { lut[m] = 1.0f; continue; }
        float J = (float)m / (float)SKETCH_S;
        float id = 1.0f + (1.0f / (float)KMER_K) * logf(2.0f * J / (1.0f + J));
        lut[m] = id < 0.0f ? 0.0f : id;
    }

    // Precompute each unit's bin so binning is a single lookup.
    int *bin = malloc(sizeof(int) * N);
    for (int i = 0; i < N; i++) bin[i] = (int)((long long)i * dim / N);

    int T = num_threads();
    if (T > MATRIX_MAX_THREADS) T = MATRIX_MAX_THREADS;
    if (T > N) T = N;
    if (T < 1) T = 1;
    if (!g_quiet) printf("  all-vs-all on %d threads ...\n", T);

    size_t cells = (size_t)dim * dim;
    double  **psum = malloc(sizeof(double *)  * T);
    uint32_t **pcnt = malloc(sizeof(uint32_t *) * T);
    pthread_t *th = malloc(sizeof(pthread_t) * T);
    MatArg    *args = malloc(sizeof(MatArg) * T);

    for (int t = 0; t < T; t++) {
        psum[t] = calloc(cells, sizeof(double));
        pcnt[t] = calloc(cells, sizeof(uint32_t));
        args[t] = (MatArg){ sig, bin, lut, N, dim, t, T, psum[t], pcnt[t] };
        pthread_create(&th[t], NULL, mat_worker, &args[t]);
    }
    for (int t = 0; t < T; t++) pthread_join(th[t], NULL);

    // Reduce private buffers into thread 0's, then average.
    double   *sum = psum[0];
    uint32_t *cnt = pcnt[0];
    for (int t = 1; t < T; t++) {
        for (size_t k = 0; k < cells; k++) { sum[k] += psum[t][k]; cnt[k] += pcnt[t][k]; }
        free(psum[t]); free(pcnt[t]);
    }

    float *mat = malloc(cells * sizeof(float));
    for (size_t k = 0; k < cells; k++)
        mat[k] = cnt[k] ? (float)(sum[k] / cnt[k]) : 0.0f;

    free(psum[0]); free(pcnt[0]);
    free(psum); free(pcnt); free(th); free(args); free(bin);
    return mat;
}

// ----------------------------------------------------------------------------
// Colormap (turbo-ish; low identity -> cool, high -> warm)
// ----------------------------------------------------------------------------

static Color colormap(float t) {
    if (t < 0.0f) t = 0.0f; else if (t > 1.0f) t = 1.0f;
    // Piecewise blue -> cyan -> green -> yellow -> red.
    float r, g, b;
    if (t < 0.25f)      { float u = t / 0.25f;        r = 0;            g = u;            b = 1; }
    else if (t < 0.50f) { float u = (t - 0.25f) / 0.25f; r = 0;         g = 1;            b = 1 - u; }
    else if (t < 0.75f) { float u = (t - 0.50f) / 0.25f; r = u;         g = 1;            b = 0; }
    else                { float u = (t - 0.75f) / 0.25f; r = 1;         g = 1 - u;        b = 0; }
    return (Color){ (unsigned char)(r * 255), (unsigned char)(g * 255), (unsigned char)(b * 255), 255 };
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file.fasta>\n", argv[0]);
        return 1;
    }

    printf("Loading %s ...\n", argv[1]);
    Fasta fa = load_fasta(argv[1]);
    printf("  %d repeat units\n", fa.count);
    if (fa.count == 0) { fprintf(stderr, "No sequences found.\n"); return 1; }

    // Headless benchmark mode: repeat the full compute pipeline and report timing.
    if (argc >= 4 && strcmp(argv[2], "--bench") == 0) {
        int iters = atoi(argv[3]);
        g_quiet = 1;
        printf("Benchmarking %d iterations (sketches + all-vs-all matrix)...\n", iters);
        double t0 = now_ms();
        int dim = 0;
        for (int it = 0; it < iters; it++) {
            uint32_t *sig = build_signatures(&fa);
            float *mat = build_matrix(&fa, sig, &dim);
            free(sig);
            free(mat);
        }
        double total = now_ms() - t0;
        printf("total: %.3f s   per iter: %.4f s   (%d units, %dx%d matrix)\n",
               total / 1000.0, total / 1000.0 / iters, fa.count, dim, dim);
        return 0;
    }

    printf("Building MinHash sketches (k=%d, s=%d) ...\n", KMER_K, SKETCH_S);
    uint32_t *sig = build_signatures(&fa);

    int dim = 0;
    float *mat = build_matrix(&fa, sig, &dim);

    // Normalize color scale to the off-diagonal identity range so structure pops.
    float lo = 1.0f, hi = 0.0f;
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            if (i != j) {
                float v = mat[(size_t)i * dim + j];
                if (v < lo) lo = v;
                if (v > hi) hi = v;
            }
    if (hi <= lo) hi = lo + 1e-4f;
    printf("Identity range (off-diagonal): %.3f - %.3f\n", lo, hi);

    // Bake the matrix into an image / texture.
    Image img = GenImageColor(dim, dim, BLACK);
    Color *px = (Color *)img.data;
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) {
            float v = mat[(size_t)i * dim + j];
            float t = (v - lo) / (hi - lo);
            px[(size_t)i * dim + j] = colormap(t);
        }

    const int W = 1000, H = 1000;
    InitWindow(W, H, "censim - stained glass (self-identity)");
    SetTargetFPS(60);
    Texture2D tex = LoadTextureFromImage(img);
    SetTextureFilter(tex, TEXTURE_FILTER_POINT);
    UnloadImage(img);

    // View state: scale = screen px per texel, offset = screen pos of texel (0,0).
    float fit = (float)(H - 120) / dim;
    float scale = fit;
    Vector2 offset = { (W - dim * scale) * 0.5f, 90.0f };

    while (!WindowShouldClose()) {
        // Zoom about the cursor.
        float wheel = GetMouseWheelMove();
        if (wheel != 0) {
            Vector2 m = GetMousePosition();
            float ns = scale * (wheel > 0 ? 1.15f : 1.0f / 1.15f);
            if (ns < fit * 0.5f) ns = fit * 0.5f;
            if (ns > 40.0f) ns = 40.0f;
            // Keep the texel under the cursor fixed.
            offset.x = m.x - (m.x - offset.x) * (ns / scale);
            offset.y = m.y - (m.y - offset.y) * (ns / scale);
            scale = ns;
        }
        // Pan.
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 d = GetMouseDelta();
            offset.x += d.x; offset.y += d.y;
        }
        if (IsKeyPressed(KEY_R)) {
            scale = fit;
            offset = (Vector2){ (W - dim * scale) * 0.5f, 90.0f };
        }

        BeginDrawing();
        ClearBackground((Color){ 18, 18, 22, 255 });
        DrawTextureEx(tex, offset, 0.0f, scale, WHITE);

        DrawText(TextFormat("%s   |   %d units   |   matrix %dx%d%s",
                            GetFileName(argv[1]), fa.count, dim, dim,
                            fa.count > dim ? " (binned)" : ""),
                 20, 20, 20, RAYWHITE);
        DrawText(TextFormat("identity %.3f - %.3f    wheel: zoom   drag: pan   R: reset",
                            lo, hi),
                 20, 48, 16, (Color){ 170, 170, 180, 255 });

        // Colorbar.
        int cbx = W - 40, cby = 90, cbw = 16, cbh = H - 180;
        for (int y = 0; y < cbh; y++) {
            float t = 1.0f - (float)y / cbh;
            DrawRectangle(cbx, cby + y, cbw, 1, colormap(t));
        }
        DrawText(TextFormat("%.2f", hi), cbx - 34, cby - 4, 12, RAYWHITE);
        DrawText(TextFormat("%.2f", lo), cbx - 34, cby + cbh - 8, 12, RAYWHITE);

        EndDrawing();
    }

    UnloadTexture(tex);
    CloseWindow();

    free(mat); free(sig);
    for (int i = 0; i < fa.count; i++) free(fa.seqs[i]);
    free(fa.seqs);
    return 0;
}

#include "stained_glass.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

// ----------------------------------------------------------------------------
// k-mer sketching
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
        default:            return -1;   // N / ambiguous breaks the k-mer window
    }
}

static void sketch_unit(const char *seq, const uint32_t *seeds, uint32_t *s) {
    for (int h = 0; h < SG_SKETCH; h++) s[h] = 0xFFFFFFFFu;
    uint32_t code = 0;
    int valid = 0;
    for (const char *p = seq; *p; p++) {
        int b = base2bit(*p);
        if (b < 0) { valid = 0; code = 0; continue; }
        code = (code << 2) | (uint32_t)b;   // k=16 => full 32-bit window, no mask
        if (++valid >= SG_KMER_K) {
            for (int h = 0; h < SG_SKETCH; h++) {
                uint32_t hv = fmix32(code ^ seeds[h]);
                if (hv < s[h]) s[h] = hv;
            }
        }
    }
}

// Low identity -> cool, high -> warm (blue->cyan->green->yellow->red).
static Color sg_colormap(float t) {
    if (t < 0.0f) t = 0.0f; else if (t > 1.0f) t = 1.0f;
    float r, g, b;
    if (t < 0.25f)      { float u = t / 0.25f;            r = 0; g = u;     b = 1; }
    else if (t < 0.50f) { float u = (t - 0.25f) / 0.25f;  r = 0; g = 1;     b = 1 - u; }
    else if (t < 0.75f) { float u = (t - 0.50f) / 0.25f;  r = u; g = 1;     b = 0; }
    else                { float u = (t - 0.75f) / 0.25f;  r = 1; g = 1 - u; b = 0; }
    return (Color){ (unsigned char)(r * 255), (unsigned char)(g * 255), (unsigned char)(b * 255), 255 };
}

// ----------------------------------------------------------------------------
// Worker compute: snapshot (sg->snap, n units) -> sg->back + res_lo/res_hi.
// Runs on the background thread only; touches no GPU / raylib state.
// ----------------------------------------------------------------------------

static void sg_compute_core(StainedGlass *sg, int n) {
    if (n < 2) return;

    uint32_t seeds[SG_SKETCH];
    for (int h = 0; h < SG_SKETCH; h++) seeds[h] = fmix32((uint32_t)h * 0x9E3779B1u + 1u);
    for (int a = 0; a < n; a++)
        sketch_unit(sg->snap + (size_t)a * SG_SEQ_STRIDE, seeds, sg->sig + (size_t)a * SG_SKETCH);

    // Identity depends only on the integer match count (0..SG_SKETCH).
    float lut[SG_SKETCH + 1];
    for (int m = 0; m <= SG_SKETCH; m++) {
        if (m == 0)         { lut[m] = 0.0f; continue; }
        if (m == SG_SKETCH) { lut[m] = 1.0f; continue; }
        float J = (float)m / (float)SG_SKETCH;
        float id = 1.0f + (1.0f / (float)SG_KMER_K) * logf(2.0f * J / (1.0f + J));
        lut[m] = id < 0.0f ? 0.0f : id;
    }

    int dim = n < SG_DIM ? n : SG_DIM;
    size_t cells = (size_t)dim * dim;
    memset(sg->avg, 0, sizeof(float) * cells);
    memset(sg->cnt, 0, sizeof(uint32_t) * cells);
    for (int i = 0; i < n; i++) sg->bin[i] = (int)((long long)i * dim / n);

    for (int i = 0; i < n; i++) {
        const uint32_t *si = sg->sig + (size_t)i * SG_SKETCH;
        int bi = sg->bin[i];
        for (int j = i; j < n; j++) {
            const uint32_t *sj = sg->sig + (size_t)j * SG_SKETCH;
            int matches = 0;
            for (int h = 0; h < SG_SKETCH; h++) matches += (si[h] == sj[h]);
            float id = lut[matches];
            int bj = sg->bin[j];
            size_t x = (size_t)bi * dim + bj;
            size_t y = (size_t)bj * dim + bi;
            sg->avg[x] += id; sg->cnt[x]++;
            if (x != y) { sg->avg[y] += id; sg->cnt[y]++; }
        }
    }

    // Average each bin; find the off-diagonal identity range.
    float cur_lo = 1.0f, cur_hi = 0.0f;
    for (int r = 0; r < dim; r++)
        for (int c = 0; c < dim; c++) {
            size_t k = (size_t)r * dim + c;
            float v = sg->cnt[k] ? sg->avg[k] / sg->cnt[k] : 0.0f;
            sg->avg[k] = v;
            if (r != c) { if (v < cur_lo) cur_lo = v; if (v > cur_hi) cur_hi = v; }
        }
    if (cur_hi <= cur_lo) cur_hi = cur_lo + 1e-4f;

    // Smooth the colour range so the panel doesn't flicker between refreshes.
    if (!sg->ema_init) { sg->ema_lo = cur_lo; sg->ema_hi = cur_hi; sg->ema_init = true; }
    else { sg->ema_lo += 0.15f * (cur_lo - sg->ema_lo); sg->ema_hi += 0.15f * (cur_hi - sg->ema_hi); }
    float lo = sg->ema_lo, hi = sg->ema_hi;
    if (hi <= lo) hi = lo + 1e-4f;
    float inv = 1.0f / (hi - lo);

    // Fill the SG_DIM back buffer, upscaling (nearest) when dim < SG_DIM.
    for (int py = 0; py < SG_DIM; py++) {
        int my = (dim == SG_DIM) ? py : (int)((long long)py * dim / SG_DIM);
        for (int px = 0; px < SG_DIM; px++) {
            int mx = (dim == SG_DIM) ? px : (int)((long long)px * dim / SG_DIM);
            float v = sg->avg[(size_t)my * dim + mx];
            sg->back[(size_t)py * SG_DIM + px] = sg_colormap((v - lo) * inv);
        }
    }
    sg->res_lo = lo;
    sg->res_hi = hi;
}

static void *sg_worker(void *arg) {
    StainedGlass *sg = arg;
    pthread_mutex_lock(&sg->mtx);
    for (;;) {
        while (sg->state != SG_REQUESTED && !sg->quit)
            pthread_cond_wait(&sg->cv, &sg->mtx);
        if (sg->quit) { pthread_mutex_unlock(&sg->mtx); return NULL; }

        sg->state = SG_COMPUTING;
        int n = sg->snap_n;
        pthread_mutex_unlock(&sg->mtx);   // compute without holding the lock

        sg_compute_core(sg, n);

        pthread_mutex_lock(&sg->mtx);
        sg->state = SG_DONE;
    }
}

// ----------------------------------------------------------------------------

void sg_init(StainedGlass *sg) {
    memset(sg, 0, sizeof(*sg));
    Image img = GenImageColor(SG_DIM, SG_DIM, (Color){ 18, 18, 22, 255 });
    sg->tex = LoadTextureFromImage(img);
    SetTextureFilter(sg->tex, TEXTURE_FILTER_BILINEAR);
    UnloadImage(img);

    sg->back = malloc(sizeof(Color)    * SG_DIM * SG_DIM);
    sg->avg  = malloc(sizeof(float)    * SG_DIM * SG_DIM);
    sg->cnt  = malloc(sizeof(uint32_t) * SG_DIM * SG_DIM);
    sg->sig  = malloc(sizeof(uint32_t) * SG_MAX_UNITS * SG_SKETCH);
    sg->bin  = malloc(sizeof(int)      * SG_MAX_UNITS);
    sg->snap = malloc((size_t)SG_MAX_UNITS * SG_SEQ_STRIDE);

    sg->last_submit = -1.0e9;
    sg->last_gen = -1;
    sg->last_units = -1;
    sg->lo = 0.0f; sg->hi = 1.0f;
    sg->state = SG_IDLE;

    pthread_mutex_init(&sg->mtx, NULL);
    pthread_cond_init(&sg->cv, NULL);
    pthread_create(&sg->thread, NULL, sg_worker, sg);
}

void sg_free(StainedGlass *sg) {
    pthread_mutex_lock(&sg->mtx);
    sg->quit = true;
    pthread_cond_signal(&sg->cv);
    pthread_mutex_unlock(&sg->mtx);
    pthread_join(sg->thread, NULL);
    pthread_mutex_destroy(&sg->mtx);
    pthread_cond_destroy(&sg->cv);

    UnloadTexture(sg->tex);
    free(sg->back); free(sg->avg); free(sg->cnt); free(sg->sig); free(sg->bin); free(sg->snap);
    memset(sg, 0, sizeof(*sg));
}

void sg_update_draw(StainedGlass *sg, const Simulation *sim, Rectangle box, double interval_s) {
    double now = GetTime();

    pthread_mutex_lock(&sg->mtx);

    // 1) Upload a finished result (GPU work stays on the main thread).
    if (sg->state == SG_DONE) {
        UpdateTexture(sg->tex, sg->back);
        sg->lo = sg->res_lo;
        sg->hi = sg->res_hi;
        sg->computed = true;
        sg->state = SG_IDLE;
    }

    // 2) Submit new work if the worker is idle, throttled, and the array changed.
    if (sg->state == SG_IDLE) {
        bool changed = (sim->stats.generation != sg->last_gen) ||
                       (sim->array.num_units  != sg->last_units);
        if (!sg->computed || (now - sg->last_submit >= interval_s && changed)) {
            int N = sim->array.num_units;
            if (N >= 2) {
                int n = N < SG_MAX_UNITS ? N : SG_MAX_UNITS;
                for (int a = 0; a < n; a++) {
                    int idx = (int)((long long)a * N / n);
                    const char *s = sim->array.units[idx];
                    char *dst = sg->snap + (size_t)a * SG_SEQ_STRIDE;
                    size_t L = 0;
                    while (s[L] && L < SG_SEQ_STRIDE - 1) { dst[L] = s[L]; L++; }
                    dst[L] = '\0';
                }
                sg->snap_n     = n;
                sg->last_submit = now;
                sg->last_gen    = sim->stats.generation;
                sg->last_units  = N;
                sg->state = SG_REQUESTED;
                pthread_cond_signal(&sg->cv);
            }
        }
    }
    pthread_mutex_unlock(&sg->mtx);

    // 3) Draw (main-owned texture + last uploaded range).
    DrawRectangleRec(box, (Color){ 18, 18, 22, 235 });
    DrawTexturePro(sg->tex, (Rectangle){ 0, 0, SG_DIM, SG_DIM }, box,
                   (Vector2){ 0, 0 }, 0.0f, WHITE);
    DrawRectangleLinesEx(box, 1.0f, (Color){ 0, 180, 0, 200 });

    DrawText("SELF-IDENTITY", (int)box.x + 6, (int)box.y + 5, 14, (Color){ 0, 200, 0, 220 });
    const char *sub = TextFormat("%d units  id %.2f-%.2f",
                                 sim->array.num_units, sg->lo, sg->hi);
    DrawText(sub, (int)box.x + 6, (int)(box.y + box.height) - 18, 12,
             (Color){ 200, 200, 210, 220 });
}

#include "reference.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void reference_init(Reference *r) {
    memset(r, 0, sizeof(*r));
}

int reference_load(Reference *r, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    reference_init(r);

    int ver = 0;
    if (fscanf(f, " CENSIM_REF %d", &ver) != 1) { fclose(f); return 0; }
    fscanf(f, " arrays %d", &r->narrays);

    char name[24];
    float mn, mx;
    int lg, nb;
    int cap = (int)(sizeof(r->metrics) / sizeof(r->metrics[0]));
    while (r->count < cap &&
           fscanf(f, " metric %23s %f %f %d %d", name, &mn, &mx, &lg, &nb) == 5) {
        if (nb <= 0 || nb > 100000) break;
        RefMetric *m = &r->metrics[r->count];
        strncpy(m->name, name, sizeof(m->name) - 1);
        m->name[sizeof(m->name) - 1] = '\0';
        hist_init(&m->h, mn, mx, nb, lg);
        long sum = 0;
        for (int i = 0; i < nb; i++) {
            long c = 0;
            if (fscanf(f, " %ld", &c) != 1) c = 0;
            m->h.counts[i] = c;
            sum += c;
        }
        m->h.total = sum;
        r->count++;
    }
    fclose(f);
    r->loaded = (r->count > 0);
    return r->loaded;
}

int reference_load_run(Reference *r, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    reference_init(r);
    int cap = (int)(sizeof(r->metrics) / sizeof(r->metrics[0]));

    char line[1024];
    // Pass 1: the "# <metric>: min=.. max=.. nbins=.. log_scale=.." headers
    // define each ghost histogram (same lines import_histograms reads).
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') {
            char name[24]; float mn, mx; int lg, nb; long tot, uf, of;
            if (sscanf(line, "# %23[^:]: min=%f max=%f nbins=%d log_scale=%d total=%ld underflow=%ld overflow=%ld",
                       name, &mn, &mx, &nb, &lg, &tot, &uf, &of) == 8) {
                if (r->count < cap && nb > 0 && nb <= 100000) {
                    RefMetric *m = &r->metrics[r->count++];
                    strncpy(m->name, name, sizeof(m->name) - 1);
                    m->name[sizeof(m->name) - 1] = '\0';
                    hist_init(&m->h, mn, mx, nb, lg);
                    m->h.total = tot; m->h.underflow = uf; m->h.overflow = of;
                }
            }
        } else if (strncmp(line, "metric,", 7) == 0) {
            break;  // column header -> data rows follow
        }
    }
    // Pass 2: "metric,bin,bin_lo,bin_center,bin_hi,count" rows.
    while (fgets(line, sizeof(line), f)) {
        char *c1 = strchr(line, ',');
        char *clast = strrchr(line, ',');
        if (!c1 || !clast || clast == c1) continue;
        int len = (int)(c1 - line); if (len > 23) len = 23;
        char metric[24]; strncpy(metric, line, len); metric[len] = '\0';
        int bin = atoi(c1 + 1);
        long count = atol(clast + 1);
        for (int i = 0; i < r->count; i++) {
            if (strcmp(r->metrics[i].name, metric) == 0) {
                if (bin >= 0 && bin < r->metrics[i].h.nbins) r->metrics[i].h.counts[bin] = count;
                break;
            }
        }
    }
    fclose(f);
    r->loaded = (r->count > 0);
    return r->loaded;
}

const Histogram *reference_get(const Reference *r, const char *name) {
    if (!r->loaded) return NULL;
    for (int i = 0; i < r->count; i++)
        if (strcmp(r->metrics[i].name, name) == 0) return &r->metrics[i].h;
    return NULL;
}

void reference_free(Reference *r) {
    for (int i = 0; i < r->count; i++) hist_free(&r->metrics[i].h);
    reference_init(r);
}

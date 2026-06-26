#include "reference.h"
#include <stdio.h>
#include <string.h>

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

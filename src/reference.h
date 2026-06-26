#ifndef REFERENCE_H
#define REFERENCE_H

#include "hist.h"

// Real-data reference distributions loaded from a reference.dat produced by
// scripts/build_reference.py (aggregated TRASH HOR tables). Each metric is a
// binned histogram whose range/scale matches the dashboard's, so it can be drawn
// as an overlay density curve. Look up by metric name (e.g. "block_gap").

typedef struct {
    char      name[24];
    Histogram h;
} RefMetric;

typedef struct {
    int       loaded;
    int       narrays;        // how many real arrays were aggregated
    int       count;
    RefMetric metrics[8];
} Reference;

void reference_init(Reference *r);
// Load reference.dat; returns 1 on success, 0 if missing/invalid.
int  reference_load(Reference *r, const char *path);
// Matching metric histogram, or NULL.
const Histogram *reference_get(const Reference *r, const char *name);
void reference_free(Reference *r);

#endif // REFERENCE_H

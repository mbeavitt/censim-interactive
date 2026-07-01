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
// Load reference.dat (CENSIM_REF format); returns 1 on success, 0 otherwise.
int  reference_load(Reference *r, const char *path);
// Load any sweep/run histogram CSV (the format export_histograms writes) as the
// reference/ghost; returns 1 on success. Lets any past run act as the overlay.
int  reference_load_run(Reference *r, const char *path);
// Matching metric histogram, or NULL.
const Histogram *reference_get(const Reference *r, const char *name);
void reference_free(Reference *r);

#endif // REFERENCE_H

# TODO

## Real reference data (deferred)

**Use real A. thaliana data for the dashboard "ghost" overlays.**

Currently the dashboard target curves are *estimated from the R1 report's summary
statistics* (means/variances/percentiles), not the real empirical distributions.
Swap these for the actual binned-histogram caches when available.

- Source: `caches.tar.gz` on the rsync.net box (`ssh rsync`). The restricted shell
  blocked `tar` over ssh; need another route to read the index (e.g. `ssh rsync help`
  to find an allowed listing tool, or pull the archive somewhere we can open it).
- Hunting for binned histogram files: names containing
  `hor`, `block`, `gap`, `similar`, `unique`, `diversit`, or `.csv`/`.tsv`/`bins`.
- Wire-in point: see `REFERENCE_*` stand-in constants in the dashboard source
  (search for `// TODO: use real data`).

### Stand-in stats currently hardcoded (from R1 report)

Groups: A. thaliana (real), Uniform Recombination (sim), Kinetochore-assoc / KARMA (sim).

| Metric | A. thaliana | Uniform | KARMA |
|---|---|---|---|
| Unique repeats / kb | mu=2.3, var=0.080 | mu=1.5, var=0.006 | mu=1.1, var=0.014 |
| HORs / kb | mu=261.4, var=97422 | mu=124.5, var=1679 | mu=192.6, var=4696 |
| HOR similarity (median/90th/99th) | 0.16 / 0.23 / 0.40 | 0.16 / 0.23 / 0.40 | excess of high values |
| Block size (median; 99th pct) | 3-4; 10 | 3-4; 85 | 3-4; 111 |
| Internal diversity (frac of HORs with score 1.0) | ~0.95 | ~0.65 | ~0.35 |
| Block gap | logarithmic; more mid/large gaps | many small gaps | fewer small, more large |
| Composite | orders of magnitude higher than both sims | baseline | baseline |

HOR detection params (match report / TRASH binary): threshold = **7** SNVs
(= floor(4% * 178), the percentage-scaled `--hor_threshold 4`), cutoff (min block len) = 3.
Validated byte-equal to `~/Code/R/TRASH_2/dep/HOR.V3.3` (-m 1); no mafft needed for
fixed-width units. See memory `hor-detection-equivalence`.

## Paper baseline (from R1 report line 46 + censim code)

| param | value |
|---|---|
| starting array | 15,000 x 178bp |
| SNP rate | Poisson 0.1 |
| INDEL rate | Poisson 0.5 |
| INDEL size | Poisson 7.6 repeats |
| dup/del | 50/50 |
| generations | 6,000,000 |
| bounding | none (unbounded drift) |
| collapse threshold | report says **2000**; censim `hpc_orig` code says 300 |
| expected collapse rate | ~15% |

DECISION: collapse threshold is user-pickable; **default 300** (per user). Trajectory
view starts at 15,000 (single view stays 10,000). The report/code disagree on the
threshold (2000 vs 300) — flagged; using 300 as default.

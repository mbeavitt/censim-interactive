#!/usr/bin/env python3
"""Build the dashboard's real-data reference from TRASH HOR tables + aligned fastas.

For each (accession, chromosome) it reads the raw HOR.V3.3 table
  HORs_method_1_aligned_<...>_Chr<N>_CEN178.fasta_t_7_c_3.csv
    columns: start_A,end_A,start_B,end_B,direction,total_variant   (1-based indices)
and the matching aligned repeat fasta
  aligned_<...>_Chr<N>_CEN178.fasta
and computes the same six metrics the dashboard plots:

  per-HOR  : block size, block gap, similarity = 1/(1+total_variant/size),
             composite = gap * similarity * size * diversity
  per-array: unique repeats / kb, HORs / kb
  (diversity = distinct repeats in block A / block size; needs the fasta)

Values are folded into fixed histograms whose ranges/scales MATCH the dashboard
(see batch.c hist_init). Run as a SLURM array (subcommand `one`, one task per
chromosome -> partial .npz), then `merge` the partials into reference.dat.

Usage:
  build_reference.py one   --hor <csv> [--fasta <fasta>] --out <partial.npz>
  build_reference.py merge --inputs <glob-or-dir> --out reference.dat
"""
import argparse, glob, os, sys
import numpy as np

REPEAT_SIZE = 178

# name -> (min, max, log_scale).  MUST match src/batch.c hist_init ranges.
SPECS = [
    ("unique_per_kb", 0.0,      6.0,   0),
    ("hors_per_kb",   1.0,      1e6,   1),
    ("block_size",    1.0,      1e6,   1),
    ("block_gap",     1.0,      1e7,   1),
    ("similarity",    0.0,      1.0,   0),
    ("composite",     1.0,      1e15,  1),
]
NBINS = 64  # reference curve resolution (independent of the sim's bin count)


def bin_indices(values, vmin, vmax, log, nbins):
    """Vectorised value -> bin index, matching the C hist_bin_index semantics:
    log axis folds <=0 and <min into bin 0; both clamp the top into the last bin."""
    v = np.asarray(values, dtype=np.float64)
    counts = np.zeros(nbins, dtype=np.int64)
    if log:
        idx = np.empty(v.shape, dtype=np.int64)
        nonpos = v <= 0
        idx[nonpos] = 0
        pos = ~nonpos
        lo, hi = np.log10(vmin), np.log10(vmax)
        lv = np.log10(np.where(pos, v, 1.0))
        frac = (lv - lo) / (hi - lo)
        j = np.floor(frac * nbins).astype(np.int64)
        j[lv < lo] = 0
        idx[pos] = j[pos]
    else:
        frac = (v - vmin) / (vmax - vmin)
        idx = np.floor(frac * nbins).astype(np.int64)
    np.clip(idx, 0, nbins - 1, out=idx)
    return np.bincount(idx, minlength=nbins).astype(np.int64)


def read_fasta_seqs(path):
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur)); cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def fasta_for_hor(hor_csv):
    d, base = os.path.split(hor_csv)
    core = base[len("HORs_method_1_"):] if base.startswith("HORs_method_1_") else base
    for suf in ("_t_7_c_3.csv", ".csv"):
        if core.endswith(suf):
            core = core[: -len(suf)]; break
    return os.path.join(d, core)  # e.g. aligned_<...>_Chr<N>_CEN178.fasta


def process_one(hor_csv, fasta):
    seqs = read_fasta_seqs(fasta)
    n = len(seqs)
    out = {name: np.zeros(NBINS, dtype=np.int64) for name, *_ in SPECS}
    if n == 0:
        return out, 0
    # integer id per distinct (aligned) repeat sequence, in file order
    _, ids = np.unique(np.array(seqs, dtype=object), return_inverse=True)
    kb = n * REPEAT_SIZE / 1000.0

    data = np.loadtxt(hor_csv, delimiter=",", skiprows=1,
                      usecols=(0, 1, 2, 3, 5), dtype=np.int64, ndmin=2)
    if data.size == 0:
        # still record the per-array metrics (0 HORs)
        out["unique_per_kb"] += bin_indices([len(np.unique(ids)) / kb], *SPECS[0][1:], NBINS)
        out["hors_per_kb"]   += bin_indices([0.0], *SPECS[1][1:], NBINS)
        return out, 1
    sA, eA, sB, eB, tv = data.T
    size = (eA - sA + 1).astype(np.float64)
    gap = np.maximum(0, sB - eA - 1).astype(np.float64)
    sim = 1.0 / (1.0 + tv / size)

    # diversity per HOR: distinct repeat ids within block A (small blocks -> cheap)
    a = sA - 1
    div = np.empty(size.shape, dtype=np.float64)
    for k in range(len(size)):
        blk = ids[a[k]: eA[k]]
        div[k] = (np.unique(blk).size) / size[k]
    composite = gap * sim * size * div

    out["block_size"] += bin_indices(size, *SPECS[2][1:], NBINS)
    out["block_gap"]  += bin_indices(gap,  *SPECS[3][1:], NBINS)
    out["similarity"] += bin_indices(sim,  *SPECS[4][1:], NBINS)
    out["composite"]  += bin_indices(composite, *SPECS[5][1:], NBINS)
    out["unique_per_kb"] += bin_indices([len(np.unique(ids)) / kb], *SPECS[0][1:], NBINS)
    out["hors_per_kb"]   += bin_indices([len(size) / kb], *SPECS[1][1:], NBINS)
    return out, 1


def cmd_one(args):
    fasta = args.fasta or fasta_for_hor(args.hor)
    if not os.path.exists(fasta):
        sys.exit(f"aligned fasta not found: {fasta}")
    out, narr = process_one(args.hor, fasta)
    np.savez(args.out, narrays=narr, **out)
    print(f"{os.path.basename(args.hor)}: arrays={narr}")


def cmd_merge(args):
    paths = []
    for p in args.inputs:
        paths += glob.glob(os.path.join(p, "*.npz")) if os.path.isdir(p) else glob.glob(p)
    if not paths:
        sys.exit("no partial .npz inputs found")
    totals = {name: np.zeros(NBINS, dtype=np.int64) for name, *_ in SPECS}
    narrays = 0
    for p in paths:
        z = np.load(p)
        narrays += int(z["narrays"])
        for name, *_ in SPECS:
            totals[name] += z[name]
    with open(args.out, "w") as f:
        f.write("CENSIM_REF 1\n")
        f.write(f"arrays {narrays}\n")
        for name, vmin, vmax, log in SPECS:
            f.write(f"metric {name} {vmin:g} {vmax:g} {log} {NBINS}\n")
            f.write(" ".join(str(int(c)) for c in totals[name]) + "\n")
    print(f"merged {len(paths)} partials ({narrays} arrays) -> {args.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    o = sub.add_parser("one"); o.add_argument("--hor", required=True)
    o.add_argument("--fasta"); o.add_argument("--out", required=True)
    o.set_defaults(func=cmd_one)
    m = sub.add_parser("merge"); m.add_argument("--inputs", nargs="+", required=True)
    m.add_argument("--out", default="reference.dat"); m.set_defaults(func=cmd_merge)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

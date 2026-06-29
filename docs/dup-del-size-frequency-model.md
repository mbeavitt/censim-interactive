# The duplication/deletion size–frequency tradeoff

Modelling assumptions behind the single **"Dup/del ratio"** control in the
trajectory and single-array views. This is a deliberately load-bearing part of
the model — it encodes a biological assumption about how centromeric satellite
arrays stay stable — so the reasoning, the maths, and the numerical caveats are
documented here for reference (and for any future methods section).

---

## 1. Biological motivation

A centromeric repeat array is a tandem run of near-identical monomers
(`REPEAT_SIZE = 178 bp` here). It evolves mainly through **unequal
crossing-over / replication slippage**, which **duplicates** or **deletes**
blocks of whole repeat units. Over evolutionary time real arrays neither shrink
to nothing nor grow without bound: they fluctuate around a characteristic size.

We want a *single knob* that lets us explore an asymmetry between duplications
and deletions — "are duplications large and rare, or small and frequent?" —
**without** that knob secretly driving the array to extinction or to runaway
growth. In other words, the asymmetry parameter should be **orthogonal to the
mean array size**.

The key biological assumption we are encoding:

> **Assumption B0.** The mean array length is approximately stationary on the
> timescales of interest; whatever sets the dup/del *asymmetry* does not also set
> the *trend* in array length.

Everything below is the consequence of taking B0 seriously.

---

## 2. The parameterisation

Two independent magnitude knobs and one asymmetry knob:

| Symbol | UI label | Meaning |
| --- | --- | --- |
| `λ` (`indel_size_lambda`) | **Mean size** | central (geometric-mean) event size, in repeat units |
| `μ` (`indel_rate`) | **INDEL rate** | expected indel events per generation |
| `r` (`dup_del_size_ratio`) | **Dup/del ratio** | dup:del **size** ratio (log slider, `0.001×`–`1000×`, centred at `1×`) |

From the single asymmetry parameter `r` we derive **four** quantities (shown live
under the slider):

```
  size_dup = λ · √r          freq_dup = P(dup) = 1 / (1 + r)
  size_del = λ / √r          freq_del = P(del) = r / (1 + r)
```

- The **size** split is geometric about `λ`: `r` and `1/r` are mirror images, and
  `size_dup / size_del = r` exactly.
- The **frequency** split is forced to be the *inverse* of the size ratio
  (`P(dup)/P(del) = 1/r`), i.e. `dup_bias = 1/(1+r)`. This coupling is the whole
  point — see §3.

Interpretation: as `r` rises, duplications get **larger but proportionally
rarer**; deletions get smaller but more frequent (and vice-versa as `r` falls).

---

## 3. Driftless coupling — derivation

Let the **expected net change in array length per indel event** be `ΔE`. With the
size/frequency split above, and writing `E[size]` for the *realised* mean applied
event size:

```
ΔE = P(dup)·E[size_dup]  −  P(del)·E[size_del]
```

**Assume** (for now) that the realised mean equals the nominal mean,
`E[size_dup] = λ√r` and `E[size_del] = λ/√r`. Then:

```
ΔE = (1/(1+r))·λ√r  −  (r/(1+r))·(λ/√r)
   = λ√r/(1+r)      −  λ√r/(1+r)
   = 0          for every r.
```

So the coupling `dup_bias = 1/(1+r)` makes the **first moment** of the
length-change exactly zero, for any asymmetry. That is precisely Assumption B0:
the ratio is **drift-neutral by construction**.

The required condition — that the *realised* mean applied event size equals the
nominal `λ` on each side — is not automatic. Three implementation details
violated it; see §5.

---

## 4. Why `r` is also (approximately) variance-neutral

Drift-neutral is necessary but not sufficient for a useful knob: if changing `r`
quietly changed the *variance* of the length random walk, it would still change
how fast arrays collapse (see §6). It does not, to leading order.

The per-event variance of the length change (using `E[Δ] ≈ 0`) is

```
Var[Δ] ≈ P(dup)·E[size_dup²] + P(del)·E[size_del²].
```

Taking just the squared-mean part `E[size²] ≈ (mean)²`:

```
   P(dup)·(λ√r)²  +  P(del)·(λ/√r)²
 = (1/(1+r))·λ²r  +  (r/(1+r))·λ²/r
 = λ²r/(1+r)      +  λ²/(1+r)
 = λ² · (r + 1)/(1 + r)
 = λ².
```

The leading term is **`λ²`, independent of `r`**. So:

> **The asymmetry `r` is both drift-neutral (mean) and, to leading order,
> variance-neutral (spread). The random-walk dynamics of the array length — and
> hence the collapse timescale — are set by the mean event size `λ` and the indel
> rate `μ`, not by the asymmetry `r`.**

This is exactly what we want: `r` explores the biological dup/del asymmetry
*independently* of array stability. (A smaller sub-leading term, the Poisson
sampling variance `≈ P(dup)·size_dup + P(del)·size_del`, does depend weakly on
`r`, but it is dominated by `λ²` for the sizes of interest.)

This result is confirmed empirically in §7: collapse frequency rises sharply with
`λ` but is essentially flat in `r`.

---

## 5. Assumptions the simulator must satisfy (and how)

The §3 derivation needs `E[realised size] = λ` on each side. The original code
broke this in three ways. Each is a place where a naive implementation silently
violates the biological assumption B0 — worth flagging in a methods section.

### A1 — Event sizes must be unbiased (no min-1 floor)

The size sampler floored events to a minimum of 1 unit (`max(X, 1)`). For a
`Poisson(λ)` draw this inflates the realised mean to `λ + e^{-λ}` (the
probability mass at 0 is pushed up to 1). The inflation is **largest when `λ` is
small** — i.e. it acts hardest on whichever side has the *smaller* events,
breaking the balance asymmetrically and producing a net drift.

**Fix:** allow a 0-size draw and treat it as a **no-op event** (the event is
simply skipped). The realised Poisson mean is then exactly `λ`. At the default
`λ = 7.6` this changes almost nothing (`e^{-7.6} ≈ 0.05 %` of events); at small
effective sizes it is the difference between balance and runaway growth.

### A2 — Large events must be sampled correctly (Poisson underflow)

Poisson draws used Knuth's product method, which compares against `e^{-λ}`. In
32-bit float, `e^{-λ}` **underflows to 0 for `λ ≳ 88`**, after which the sampler
silently caps/corrupts the size of the large rare event. This is a hard numerical
cliff: e.g. with mean size 22 it triggers as the deletion mean crosses ~104,
around `r ≈ 0.045`, exactly where the array flips from stable to runaway.

**Fix:** for `λ > 30` use the Normal approximation `Poisson(λ) ≈ N(λ, λ)`
(rounded, clamped ≥ 0). This keeps the mean at `λ` with no underflow, and is
`O(1)` instead of `O(λ)`.

### A3 — Placement must be size-unbiased (no reject-on-overrun)

Events were placed by drawing a uniform start position and **rejecting** any
event whose end ran past the array boundary. The rejection probability is
proportional to event size, so the **larger** of the two event types is dropped
more often — again an asymmetric, size-dependent drift that *no* choice of
frequency can cancel (it depends on the current array length `N`, not on `r`).

**Fix:** draw the start position from the *valid* range `[0, N·REPEAT_SIZE −
span)` so the event always fits and is always applied. An event larger than the
entire array still cannot be placed — but that only happens near collapse, where
it harmlessly acts as a floor.

> **Caveat for a methods section.** A1/A2/A3 make the coupling drift-free for the
> default **Poisson** size distribution. The optional **geometric** and
> **power-law** size distributions retain a minimum size of 1 by construction, so
> their realised means differ slightly from `λ`; the coupling is then only
> approximately drift-free for those models.

---

## 6. What the coupling does *not* do (the fundamental limit)

The coupling nulls the **mean** drift. It does **not** make the array length
constant, because length is a **random walk** with non-zero variance (§4). Under
free drift (`elasticity = 0`, the paper regime) a zero-mean random walk **will**
eventually hit the collapse threshold — this is the intended drift-to-collapse
mechanism, not a bug. The expected time to collapse scales with the variance, set
by `λ` (and `μ`), **independently of `r`**.

> **To hold the array size genuinely constant you need a restoring force
> (elastic bounding), not just the frequency coupling.** The coupling guarantees
> the asymmetry knob does not *bias* the walk; elastic bounding is what removes
> the walk's diffusion. The two compose cleanly: with a drift-free base, elastic
> bounding cooperates with it instead of fighting a biased base.

So the honest statement of the model is:

- `r` (asymmetry): drift-neutral and ~variance-neutral. Safe to sweep.
- `λ`, `μ` (magnitude/rate): set the diffusion and hence the collapse timescale.
- `elasticity` (bounding): the only thing that actually pins the mean size.

---

## 7. Empirical validation

Free drift (`elasticity = 0`), start size 15 000 units, collapse threshold 300,
2 000 000 generations, mean over independent seeds. **After** the A1–A3 fixes:

Mean size `λ = 7.6` (24 seeds): mean final size of survivors stays in
`14 500 – 19 400` units across the whole slider range, with collapses `0–5 / 24`
and **no trend in `r`**:

| `r`   | 0.001 | 0.045 | 0.2 | 1 | 5 | 21 | 1000 |
| ----- | ----- | ----- | --- | - | - | -- | ---- |
| final | 14575 | 15543 | 14915 | 16203 | 19367 | 15940 | 14834 |
| coll. | 2 | 1 | 2 | 0 | 4 | 5 | 2 |

Larger mean sizes stay **bounded** (tens of thousands, never the pre-fix
millions), and collapse frequency rises with `λ` — *not* with `r`. Critically, at
`r = 1` (a perfectly symmetric dup/del process) the collapse rate is already as
high as at the extremes, confirming the collapses are the free-drift random walk
(§6), governed by `λ`, and **not** an artefact of the asymmetry:

| `λ` | `r=0.044` | `r=0.2` | `r=1` | `r=21` | `r=50` |
| --- | --------- | ------- | ----- | ------ | ------ |
| 22  | 11/16 | 5/16 | 7/16 | 11/16 | 12/16 |
| 50  | 8/16  | 13/16 | 12/16 | — | — |

**Before** the fixes, the same sweep produced runaway growth to millions of units
at small `r` (and at large `λ`), and near-total instant collapse at large `r` —
i.e. the asymmetry knob was dominating the array dynamics, violating B0.

---

## 8. One-paragraph summary for a methods section

> Duplication and deletion events are coupled through a single asymmetry
> parameter `r`: mean event sizes are split geometrically about a central size
> `λ` as `λ√r` (duplications) and `λ/√r` (deletions), and their relative
> frequencies are set to the inverse ratio, `P(dup) = 1/(1+r)`. This makes the
> expected per-event change in array length identically zero for all `r`, and
> leaves the per-event variance approximately `λ²` independent of `r`, so the
> asymmetry is explored without biasing the mean array length or its random-walk
> diffusion. Achieving this in practice requires unbiased event-size sampling
> (no minimum-size floor; a Normal approximation for large Poisson means to avoid
> single-precision underflow) and size-unbiased placement (sampling the event
> start within the valid range rather than rejecting boundary overruns).
> Under free drift the array length remains a zero-mean random walk that
> eventually collapses on a timescale set by `λ` and the indel rate; a stationary
> mean size additionally requires an explicit restoring (elastic-bounding) term.

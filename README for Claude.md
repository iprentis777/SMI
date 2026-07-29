# README for Claude

Handoff for the next agent working on SMI fODF tractography through edema.

This **replaces** the earlier version of this file, which was written when the
modulation work stalled. Most of that content survives below, but **three of its
conclusions have since been measured to be wrong** and are corrected in place.
The old text is in git history (`409e3ca`) if you need it.

**Read sections 1 and 2 before proposing anything.** Section 2 exists because
the previous handoff sent a reader after the wrong cause, and that cost a day.

---

## 1. Status

The framework is essentially built. What remains is documentation and judgement
calls, not construction.

### On master

| feature | option | default | report |
|---|---|---|---|
| Regularized deconvolution | `options.fODF_regularization` | off | `REPORT_fODF_regularization_sweep.md` |
| Anisotropy modulation | `options.fODF_modulation` | off | `REPORT_fODF_modulation.md` |
| **Post hoc outlier cap** | `options.fODF_outlier` | off | `REPORT_fODF_outlier_cap.md` |

All three are opt-in, and `out.plm`, `out.pl`, `out.kernel` are bit identical
whether any of them is on or off. That invariant is tested, not asserted — keep
it. Patches `0001`-`0005` record each change and apply with `git am`.

### On a branch, not merged

`claude/freewater-simulations` (`1be06ea`, 3 commits) — the free-water
comparison package, `REPORT_fODF_freewater.md`,
`REPORT_fODF_compartment_split.md`, two figures. **No PR yet.** It was split out
of PR #5 deliberately so the cap could merge alone. Every measured number quoted
in `REPORT_fODF_outlier_cap.md` sections 4-5 comes from this branch, so it
should land or those numbers lose their provenance.

### Never run on real data

Everything. No result in any report has touched a patient scan. Two driver
scripts live only on the user's machine and are sent as files:
`run_smi_hcp.m` (single subject) and `run_smi_batch_mod.m` (5-subject edema
batch). Both were updated this session for the cap; `run_smi_batch_mod.m` runs
`Lmax = [0 2 2 2 4 4]`, i.e. **max(Lmax) = 4**, which changes the ceiling
arithmetic — see section 6.

---

## 2. Corrections to the previous handoff

**Do not trust the old §7.6 or `REPORT_fODF_modulation.md` §3 on Tikhonov.**

1. **`lambda_tikhonov` does not damage the high `l` bands, and does not do
   much of anything.** Swept 0 to 0.8 at three noise levels: the 45 degree
   crossing is never resolved at any value, the 60 degree error moves without
   trend, the CSF peak changes by ~3%. `REPORT_fODF_modulation.md` §3
   attributed the `pl4` loss to it. That is wrong. Measured noise-free at
   `lambda_tikhonov = 0`, SMI still returns `p6` at **28%** of truth and `p4`
   at **62%** — the loss is real but the cause is the non-negativity constraint
   plus error in the estimated kernel, whose `K_l` at high `l` is tiny and very
   sensitive.

2. **The regularizer that is load-bearing against blow-ups is non-negativity,
   not Tikhonov.** The two had never been separated. At
   `lambda_tikhonov = 0` with `lambda_nonneg = 10`, SNR 15, the CSF peak is
   **0.39** — not the 1e13 the old §7.6 recorded for "no regularization",
   because non-negativity was on in that run too. Turn non-negativity off and
   at SNR 15 the CSF peak reaches 1.87 with WM/CSF contrast **0.61**, i.e.
   white matter dimmer than CSF.

3. **`lambda_nonneg = 10` is what closes the 45 degree crossing, and it is not
   on the Pareto front.** With the constraint off, SMI resolves 45 degrees
   **40/40** at **1.69 deg**, exactly the band-limited ground truth's own error.
   From 3 upward it never resolves. And `lambda_nonneg = 3` matches or beats the
   shipped 10 on WM/CSF contrast and 60 degree angular error at both SNR 30 and
   15, giving up nothing since both lose 45 degrees anyway.
   **The default was NOT changed.** It wants a cross-check against
   `example_fODF_regularization_sweep.m`'s own reconstruction-error metric
   across `lambda_nonneg` in {1, 3, 10} first. That is the single most concrete
   open task in the repo and is about an hour of work.

4. **§7.5's rejection of tissue-fraction weights is bounded, not universal.**
   It was measured against an edema class with `f` collapsed to 0.30. When edema
   is modelled instead as *redistribution* — `f` held at 0.60 while
   extra-axonal water becomes free water — `f` is recovered flat across the
   whole trajectory and an `f` weight suppresses CSF completely while leaving
   edema at or above healthy amplitude. So: **a tissue-fraction weight fails on
   axonal loss and is safe under water redistribution.** Which regime real
   peritumoral tissue is in decides whether `f` is usable, and that is
   measurable in a known ROI. This is now the most promising untested weight.

---

## 3. Conventions that must not be got wrong

Pinned analytically and verified numerically (round-trip `max|err| = 1.2e-15`).

`out.plm` is in the **normalized** convention `p_00 = 1`, covering `l = 2..Lmax`
only (the `l=0` term is not stored).

```
SH coefficients :  f_lm = plm .* sqrt((2l+1)/(4*pi))
l = 0 term      :  f_00 = 1/sqrt(4*pi)        (constant in every voxel)
fODF amplitude  :  A(u) = 1/(4*pi) + sum_{l>=2,m} f_lm Y_lm(u)
Lmax from count :  Lmax = sqrt(2*Nlm + 9/4) - 3/2
```

Consequences that keep mattering:

1. **The fODF integrates to 1 in every voxel.** It carries no density
   information at all. A CSF voxel and a coherent WM voxel have equal mass.
2. **The isotropic floor is a fixed `1/(4*pi) = 0.0796`**, above MRtrix's
   default `iFOD2 -cutoff 0.05`. An unmodulated SMI fODF therefore passes the
   tractography termination test **in every voxel of the brain**, CSF included.
3. To rescale a peak to target `T`, scale the anisotropic part by
   `s = (T - 1/(4pi))/(peak - 1/(4pi))`, never `T/peak`. Both the modulation and
   the outlier cap depend on this.
4. Scaling all `l>=2` coefficients uniformly preserves peak orientation exactly.
   Clipping SH coefficients **independently** does not and can swing peak
   orientations. Never do the latter.
5. **The hard physical maximum of a band-limited fODF is
   `sum_{l even <= Lmax} (2l+1)/(4*pi)`**: 1.194 at Lmax 4, 2.228 at Lmax 6,
   3.581 at Lmax 8. Any amplitude ceiling should be compared against this and
   against what SMI actually recovers, which is much lower.

**Two different p2 exist and they are not the same quantity.** Same for p4.

- `out.kernel(:,:,:,ip2)` — fitted jointly with the kernel by polynomial
  regression on rotational invariants. `ip2 = 6` without T2 fitting, `8` with.
  `SMI.grab_kernel_pl` resolves the ambiguity from `out.shells` row 4 (TE).
- `out.pl(:,:,:,1)` — the l=2 invariant of the **deconvolved** fODF, a raw
  Euclidean norm, unclipped.

With `Lmax = [0 8 8 8]`, `RotInv_Lmax` defaults to 4, so **`out.kernel` has 7
maps, not 6**: `[f Da Depar Deperp fw p2 p4]`.

---

## 4. Environment and validation

Linux container, **GNU Octave 8/9, no MATLAB**. Working setup:

```
apt-get update && apt-get install -y octave octave-statistics octave-image
pkg load statistics; pkg load image;
```

`SMI.fit` **does run end-to-end in Octave** with three shims. They are NOT in
the repo — recreate them:

| shim | why |
|---|---|
| `round(x,n)` | MATLAB two-arg form, used in `Group_dwi_in_shells_b_beta_TE` |
| `discretize(x,edges)` | not in Octave, used in `StandardModel_MLfit_RotInvs` |
| `datetime()` | not in Octave, used only by the log writer |

Other traps, all hit at least once:

- **Octave cannot call functions defined at the end of a script**, MATLAB
  requires them there. That is why helper files are separate.
- **`SMI.vectorize` takes a different branch if any spatial dimension is a
  singleton.** Always build simulation volumes with all three dims > 1.
- **`end` inside a comma-list expansion** (`fod(p{:}, 2:end)` where `p` is a
  `num2cell` result) is miscomputed in Octave. Use explicit indices.
- **`median(...,4,'omitnan')`** is version-dependent. `SMI.neighbour_median`
  deliberately takes the median by sorting instead (NaN sorts last in both), so
  it needs no nanflag and no toolbox.
- Graphics are unusable. Stub `figure/histogram/subplot/...` to test plotting.
- `strel('disk',r,0)` — Octave requires the explicit `N`; `N=0` is an exact
  disk and is the better choice in MATLAB too.
- Python side of the comparison work needs `numpy scipy dipy cvxpy matplotlib`.
  **cvxpy is not optional** — dipy's MSMT-CSD QP fails without it.

**Parse-checking is not enough.** Runtime errors reached the user twice that
way, because the interesting code sits inside `try` blocks. Two patterns that
do work, both used this session:

1. Extract the post-fit section of the driver script, stub the NIfTI I/O, feed a
   synthetic `out`, run all flag combinations.
2. For structural edits, a **bracket-depth-aware** keyword balance check.
   A naive one is useless: `Pp(2:end-1,...)` and `addpath(..., '-end')` both
   produce false positives.

---

## 5. What is in the repo, and what each thing is for

| file | what |
|---|---|
| `SMI.m` | the toolbox. All three opt-in features live here as static methods |
| `README.md` | user-facing documentation of all three features |
| `REPORT_fODF_regularization_sweep.md` | the `lambda_nonneg = 10` measurement |
| `REPORT_fODF_modulation.md` | the anisotropy weight measurement. **§3 is wrong on Tikhonov, see section 2** |
| `REPORT_fODF_outlier_cap.md` | the cap measurement |
| `example_fODF_regularization*.m` | regularization examples and the sweep |
| `example_fODF_modulation.m` + `fODF_modulation_helpers.m` | the 7-class simulation. The helpers are the reusable forward model and are used by everything |
| `test_fODF_outlier_cap.m`, `test_SMI_outlier_cap.m` | the cap's tests. Self-contained, no data needed |
| `0001`-`0005*.patch` | one patch per change, `git am`-able |

On `claude/freewater-simulations` only: `freewater_comparison/` (the SMI vs
CSD vs MSMT-CSD package), `REPORT_fODF_freewater.md`,
`REPORT_fODF_compartment_split.md`, two figures.

---

## 6. What the three measurement campaigns found

### 6.1 SMI vs CSD vs MSMT-CSD in free water

Same synthesised DWI fed to all three, one shared evaluation sphere, one shared
peak finder, all at Lmax 6.

**SMI's fODF amplitude is blind to free water; CSD's and MSMT-CSD's scale with
the tissue fraction.** Peak ratio at `fw = 0.40` over `fw = 0`, dilution model:
SMI 0.94-1.01, CSD 0.60, MSMT 0.62. Structural, not a tuning artefact — SMI's
fODF has unit mass by construction.

The flip side, and the whole tractography problem in one line: SMI leaves CSF at
**0.32** (four times its isotropic floor, 100% above MRtrix's cutoff) where
MSMT-CSD leaves it at **0.028** and 0%. **SMI cannot dim in edema and cannot
terminate in CSF; MSMT-CSD does both.** On orientation SMI matches MSMT-CSD and
beats CSD.

### 6.2 Edema as redistribution, not dilution

Hold `f = 0.60` and convert extra-axonal water to free water. This matters more
than the free water fraction does: the same `fw = 0.40` costs CSD 13% and MSMT
23% of amplitude instead of 40%, because an AFD-style fODF tracks *fibre*
content and `f` is what did not move. See section 2, item 4 for the consequence.

Also: **MSMT-CSD's CSF volume fraction is not a usable free-water estimate under
this model** — 0.019 where the truth is 0.20. SMI's `fw` reads 0.151.

### 6.3 Response estimation

CSD and MSMT were given responses estimated from the most anisotropic voxels of
a brain containing 15% edema, rather than handed the truth.

Edema with intact axons is **not distinguishable** from healthy WM on a high-b
anisotropy criterion (b=3 anisotropy 0.677 vs 0.682), so it is enriched in the
selection — but the estimated response is unchanged to three decimals either
way. The contamination is benign because the contaminant looks like the thing
being estimated.

What is **not** benign: any estimated response is far blunter than an idealised
delta response, because it absorbs fibre dispersion. Single-shell CSD absorbs
that; **MSMT-CSD degrades 3-4x on crossings**. Note also that handing CSD an
exact delta response, as `REPORT_fODF_freewater.md` originally framed it, is
*idealised*, not generous.

---

## 7. Dead ends, with evidence

Still valid from the previous handoff. Do not re-walk these.

- **`p2product` (the shipped modulation default) destroys symmetric 3-way
  crossings** — 100% below cutoff — because `p2` is *identically zero* for three
  equal orthogonal fibres. Realistic 3-fibre voxels survive but are dimmed 4x.
- **The kernel `p4` is structurally blind to crossings.**
  `Get_uniformly_distributed_SM_prior` draws `p4 = rand * p2 * 0.9`, so with
  `p2 = 0` the regression cannot represent `p4 = 0.54`. Measured kernel
  anisotropic power for an orthogonal crossing: 0.000.
- **Raw anisotropic power fails; noise does not cancel.** `p_l` is a norm of
  noisy coefficients, positive when the truth is zero. At SNR 15 it *inverts*:
  CSF 0.207 > orthogonal crossing 0.170.
- **The signal-domain route does not escape the noise floor.** Orthogonal
  crossing 0.0447 vs CSF 0.0365. Changing domain does not create information.
- **Lmax 8 is worse than Lmax 6** on unweighted FP (52.5% vs 21.9%), and at
  Lmax 8 a genuine 3-way crossing is *dimmer than CSF noise* before any
  weighting.
- **`degenerate = 'clip'` is a bad modulation default.** In a blown-up voxel the
  raw `p` exceeds the clip so the voxel gets weight exactly 1.0, the maximum in
  the volume. `'reject'` exists. Still not changed.

The one thing that has ever satisfied every constraint at once, still unshipped:
**noise-floor-subtracted anisotropic power over `l = 2,4`**, using `sigma` and
the kernel gains already in `out`. 0.0% FP on GM+CSF at both SNRs with all six
WM classes retained. Blocked on a free safety factor nobody has justified, and
it has never touched real data.

---

## 8. Next steps, roughly in order of value per hour

1. **Run the drivers on real data and read the peak histogram.** Everything
   above is simulation. The single most informative number is
   `out.fODF_outlier.Ncap`: simulation says it should be **0** on a regularized
   fit, so a non-zero count means the regularization is not behaving on real
   data as it does in simulation, which is a bigger finding than the cap. If it
   is non-zero, look at `SMI_fODFcap_flagged.nii` — ventricles means the
   deconvolution is blowing up; brain edge means the pre-fit erosion is too
   gentle.
2. **The `lambda_nonneg` cross-check** (section 2, item 3). Concrete, bounded,
   and would justify a default change with two independent measurements.
3. **Test `f` as a modulation weight in a known edema ROI** (section 2, item 4).
   This is the question that decides whether the whole modulation approach can
   use tissue fractions after all.
4. **Verify SMI's SH basis against MRtrix's.** Still never done. Synthesise a
   delta along a known direction, write it, run `sh2peaks`, compare. Half an
   hour, and every peak orientation shipped to `tckgen` depends on it.
5. **Land `claude/freewater-simulations`**, or the measured numbers quoted in
   `REPORT_fODF_outlier_cap.md` have no reproducible source.
6. Spatial context beyond the cap — a neighbour-agreement *weight* rather than
   just an outlier test — is the remaining unexplored axis. The cap proves the
   plumbing works and that a neighbourhood statistic is inherently edema-safe.

---

## 9. Documentation and clutter debt

This is what the user asked to focus on next, and it is real.

- **`REPORT_fODF_modulation.md` §3 states a wrong cause** for the `pl4` loss.
  Fix it in place rather than leaving a corrected copy elsewhere.
- **Five reports, no index.** `README.md` links some. A reader has no map of
  which report supersedes which, or which numbers are simulation-only.
  Everything in every report is simulation-only; that is worth saying once,
  loudly, in `README.md`.
- **`0005fODFoutliercap.patch` duplicates its own diff** (~1100 of 2051 added
  lines in PR #5). That is the repo convention, so it stays, but be aware the
  patch files roughly double the apparent size of every change.
- **Removed this session:** `.DS_Store` and `SMI.asv` (a MATLAB autosave) were
  tracked. Both are now gitignored along with `log_SMI_*.txt`, which `SMI.fit`
  drops in the working directory on every run.
- **Stale branches on origin:** `add-claude-github-actions-1785176603102`,
  `claude/github-repo-markdown-changes-t49gi0`,
  `claude/fodf-modulation-edema-38eh81` are all merged and can be deleted.
  `claude/freewater-simulations` is the only unmerged one that matters.
- **CI runs `claude-review` only.** There is no MATLAB or Octave runner, so
  `test_fODF_outlier_cap.m` and `test_SMI_outlier_cap.m` are never executed by
  CI — they pass locally under Octave and that is invisible on a PR. An Octave
  job would be straightforward and the repo has never had one.

---

## 10. Working preferences observed

- Wants to be told when a suggestion is wrong, **with evidence**. Has pushed
  back correctly several times, and was right to reject tissue-type criteria
  three times before the redistribution result bounded that objection.
- **Prefers measured over assumed.** Every default in this repo that was
  guessed has turned out wrong by roughly 2x. The regularization sweep, the
  crossing finding, the ceiling, and the Tikhonov correction all exist because
  someone checked instead of assuming.
- **Nothing goes in the repo until it is proven.** Stated explicitly. Honour it.
- Asks where changes were made — point at file and line numbers.
- Wants patch files and README updates alongside every change, so the history
  can be followed without reading diffs.
- Figures matter; a manuscript is the eventual target. The
  manuscript-grade simulation is still blocked on ground-truth *density*: every
  voxel is a single kernel, so density-weighted scoring is meaningless until
  there is a second tissue type per voxel. `REPORT_fODF_modulation.md` §2
  already has 50/50 interface classes, so the machinery exists — it needs a
  continuous mixing fraction, not new code.
- Two driver scripts live only on the user's machine and arrive as file uploads.
  Peak truncation must run **before** modulation, and the outlier cap now
  supersedes it.

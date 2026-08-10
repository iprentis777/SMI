# Notebooks

Narrated, step-by-step versions of the simulations, for reading rather than
running in bulk. They are the answer to "how do I check this myself?".

| notebook | what it walks through |
|---|---|
| `smi_simulation_walkthrough.m` | the SMI arm end to end on a real HCP protocol: ground truth, kernel, forward convolution, Rician noise, `SMI.fit` at Lmax 4/6/8, peaks, MRtrix export. Conditions: single fibre, 30, 45, 60 degrees. **One SNR** |
| `smi_manuscript_60deg.m` | the manuscript cut — **60 degree crossing only**, **two kernels back to back**, **swept across SNR** — six isometric figures |
| `smi_manuscript_60deg.ipynb` | **the three-arm version.** Same experiment, healthy kernel only, with **SSST-CSD and MSMT-CSD actually running** through MRtrix3. A Jupyter notebook on the Octave kernel |
| `csd_manuscript_60deg.ipynb` | **the Python CSD arm.** The CSD half only, on a **Python** kernel, with a **dispersion-matched** response, exporting its results for the `.m` side to draw figures from |

## `csd_manuscript_60deg.ipynb` — the Python CSD arm

The division of labour the other notebooks could not have: **Python builds the
signal and drives the CSD arms; the `.m` file fits SMI and draws the figures.**
It does not fit SMI at all, which is why it runs in minutes where
`smi_manuscript_60deg.m` runs in hours.

| arm | what it is |
|---|---|
| SSST-CSD | `dwi2fod csd` on the top shell — the MRtrix3 binary |
| MSMT-CSD | `dwi2fod msmt_csd` on all four shells, 3 tissues — the binary |

### The response is dispersion matched, and that is the point

Every other arm in this package is handed the **delta** response
`r_l(b) = K_l(b) sqrt((2l+1) 4 pi)`, the response of a fibre with no dispersion.
But no fibre in this simulation is a delta: each population is a Watson at
`kappa = 16`. So this notebook builds

```
r_l(b) = K_l(b) * p_l^Watson(kappa) * sqrt((2l+1) * 4*pi)
```

which is the response of *one such population* — what a perfect `dwi2response`
run on this data would recover. `RESPONSE_MODE` switches between the two.

Measured at b = 3, normalised to `l = 0`:

| response | `r_2` | `r_4` | `r_6` |
|---|---|---|---|
| delta (exact kernel) | -0.828 | 0.441 | -0.179 |
| **dispersion matched, `kappa = 16`** | **-0.748** | **0.315** | **-0.089** |
| `dwi2response dhollander` | -0.703 | 0.342 | -0.107 |
| `dwi2response tournier` | -0.714 | 0.349 | -0.116 |
| `dwi2response fa` | -0.718 | 0.352 | -0.121 |

At `l = 2` the dispersion-matched response lands **between** the delta and the
estimators, accounting for about **69%** of the gap between them — so most of
what makes an estimated response blunter than the kernel is fibre dispersion,
which is the thing section 6.3 of the Monte Carlo report attributed it to. At
`l = 4` and `l = 6` it overshoots slightly, i.e. `kappa = 16` disperses a little
*more* than the phantom the estimators were run on. Read that comparison as
indicative only: those three were estimated on the **superseded synthetic
protocol**, not this acquisition.

`Figures/csd_response_dispersion.png` is this comparison drawn. It depends only
on the kernel and `kappa`, not on `NREP`, so it is exact rather than a smoke-run
artefact.

### Identity with the `.m` file is measured, not asserted

`dump_reference.m` writes the Octave side of the forward model into `data/`, and
`check_python_vs_octave.py` recomputes all of it in Python. Step 0 of the
notebook runs that comparison and fails loudly if it does not pass. Measured,
22 of 22 arrays agree:

| quantity | max abs err |
|---|---|
| SH basis, Lmax 8, both `CS_phase` conventions | 5.5e-15 |
| ground truth `plm` | 1.0e-15 |
| kernel invariants `K_l(b)` | 1.6e-15 |
| **noise-free signal** | **6.7e-15** |

Two rows are compared *relatively* and the reasons are worth knowing:
Octave's `textscan` and Python's `float()` round decimal strings like `2.99`
differently in the last bit (80 of 288 b values, 1.5e-16 relative), and the
Watson amplitudes are `exp(kappa) ~ 9e6`, where an absolute tolerance is a
demand for 22 significant digits.

**The noise realisations are deliberately NOT identical.** Octave's
`randn('seed',...)` legacy generator and numpy's PCG64 are different streams and
cannot be reconciled without reimplementing one inside the other. So the noisy
signal is **exported** instead: `C.signal` is what the CSD arms were actually
given, and fitting SMI to it puts all three arms on the same realisations.

### Handing the results to the `.m` side

```matlab
run('oct_path.m');
C = read_csd_export();     % arms, scores, SH coefficients, signal, config
figures_csd_arms();        % the figures, or figures_csd_arms('../Figures')
```

`test_csd_roundtrip.m` re-scores the exported SH coefficients with the **Octave**
peak finder (`helpers/fODF_peak_score.m`) and compares against what Python
printed. Measured: peak counts and spurious counts reproduce **exactly**, mean
angular error to **4.6e-10 deg**. That is what makes `C.scores` usable directly
rather than something the `.m` side has to re-derive.

### Findings, not failures

**MSMT-CSD's noise-free angular error is far above the ceiling, and a blunter
response makes it worse.** At 60 degrees, Lmax 6, `SNR = inf` — no noise at all:

| response | MSMT peak errors | SSST peak errors |
|---|---|---|
| dispersion matched | **7.80°**, 3.54° | 1.27°, 1.72° |
| delta | 2.82°, 3.54° | 1.27°, 1.72° |

SSST-CSD sits exactly on the ceiling either way. MSMT finds both lobes in both
cases — the *count* is right — but its primary lobe is displaced, and the
displacement roughly triples when the response is blunted. This is not a setup
bug: it reproduces from a third direction the report's finding that MSMT's 60
degree error is **response-limited** (6.85° estimated, 2.36° with the exact
response, section 6.4).

**The angular error uses the largest peak only, and a symmetric crossing is a
near-tie.** MSMT's two lobes differ by 7% in amplitude; which is "primary" can
flip on the last few ulps, and the reported error then jumps between the two
lobes' errors with no real change. Read it alongside the peak and spurious
counts. This is inherited deliberately — it is what `smi_manuscript_60deg.m`
scores, and changing it would make the arms incomparable.

### Its figures render, and so do the Octave ones

Both were verified in this container, which is new: `figures_csd_arms.m` writes
all three PNGs through Octave's `print` without error. **The claim elsewhere in
this repository that Octave's `print` is broken here is no longer true** for
these figures — gnuplot still ignores `camlight` and `lighting`, so glyphs are
flat shaded, but they render.

`SMOKE_TEST = True` ships as the default: one Lmax, three SNRs, `NREP = 27`, a
few minutes, 23 CHECKs. `False` gives the manuscript configuration.

## `smi_manuscript_60deg.ipynb` — the three-arm notebook

The same experiment as the `.m`, with the CSD arms wired in rather than stubbed,
and cut to the **healthy kernel only**; the edema arm waits on a synthetic
response.

| arm | what it is |
|---|---|
| SMI | `SMI.fit`, `fODF_regularization.flag_nonneg = 1`, `CS_phase = 0` |
| SSST-CSD | `dwi2fod csd` on the top shell — the MRtrix3 binary |
| MSMT-CSD | `dwi2fod msmt_csd` on all four shells, 3 tissues — the binary |

**All three responses are the exact analytic kernel response**,
`r_l(b) = K_l(b) sqrt((2l+1) 4 pi)`, so response estimation is not a confound
and a difference between arms is the deconvolution. Swapping in an estimated
response is one `RH.read_response` call in Step 6b — but see the warning there
about which acquisition a stored response was estimated on.

Three things the `.m` file's commented hook asserted that are **not true**, all
corrected in the notebook:

- **The arms do not "appear automatically" in the figures.** The scoring arrays
  and every figure were indexed by `(Lmax, SNR, condition)` with no arm
  dimension at all. There is one now.
- **The peak finder was not scale-free.** It subtracted the constant
  `1/(4*pi)`, which is the isotropic part only for an fODF with `p_00 = 1` in
  every voxel. An MRtrix FOD is unnormalised and its `l=0` term varies per
  voxel and carries apparent fibre density. The notebook subtracts **each
  voxel's own `l = 0` term**, and checks that this agrees with the constant for
  the SMI arm.
- **The response must be evaluated at MRtrix's own per-shell b**, read back with
  `mrinfo -shell_bvalues`, not at the nominal 0/1/2/3. This protocol's shells
  jitter — MRtrix reports `[0, 998.28, 1998.17, 2996.06]` s/mm² — so the
  nominal values attach each response row to a b nobody acquired at.

`SMOKE_TEST = true` is the shipped default: one Lmax, three SNRs, `NREP = 27`,
a few minutes. `false` gives the manuscript configuration.

**Its figures have been verified to render**, which is not true of anything else
in this repository. The Octave Jupyter kernel produces inline PNGs, so Figure 1's
shared radial scale is confirmed working rather than assumed. Only the `gnuplot`
toolkit is available in this container, and it ignores `camlight` and `lighting`,
so glyphs are flat-shaded here and correct in MATLAB. Each figure is wrapped in
`try`/`catch` and prints `** FIGURE FAILED **` rather than taking the run down.

`smi_manuscript_60deg.m` keeps **every knob in one Configuration block at the
top**, including the kernels, so retuning it never means opening another file.
It simulates **both tissues in one run**, back to back off the same ground
truth, protocol and noise seeds — only the kernel differs:

| preset | `[f Da Depar Deperp fw]` | extra-axonal |
|---|---|---|
| `healthy` | `[0.60 2.0 2.0 0.50 0.02]` | 0.38 |
| `edema` | `[0.10 2.4 2.7 1.15 0.35]` | 0.55 |

The `edema` kernel is shaped from a *fitted* kernel, which is why the
diffusivities are not round and why `Depar` exceeds `Da`. Two values depart from
the fit as supplied, both deliberately: **`f` is 0.10 rather than 0.05**, because
0.05 is exactly the lower bound of SMI's default training prior and a truth
sitting on the prior boundary is estimated with a one-sided bias that is
indistinguishable from a real effect; and **`fw = 0.35`**, which the fit did not
specify, modelling edema as added free water. `Deperp = 1.15` remains close to
its own prior cap of 1.2 — left as fitted, but the first place to look if the
extra-axonal diffusivities come back biased low.

Only three stateless utilities still come from `mc_config.m` — `pick_grid`,
`rotate_about`, `load_protocol_file`. They are shared rather than copied so the
fibre-axis convention and the protocol reader cannot drift from
`gen_montecarlo.m`. Changing a knob in the notebook does **not** change the
pipeline, which keeps its own values in `mc_config.m`.

`smi_simulation_walkthrough.m` draws four figures: the ground truth fibre
geometry, the kernel as a response function, the forward convolution with and
without noise, and the reconstructed fODFs at each Lmax.

`smi_manuscript_60deg.m` draws six, and every 3D panel opens in an **isometric**
view so all subfigures are readable without touching the camera:

| figure | layout |
|---|---|
| 1 | the ground truth fODF, then **each kernel's response glyph per shell**, then a **third row holding the difference between them**. All response glyphs share **one radial scale**, so size carries meaning: a smaller glyph is a genuinely smaller signal |
| 2 | the signal as surfaces on the sphere, **healthy**: b down the rows, SNR across |
| 3 | the same for **edema**, on the same scale as Figure 2 |
| 4 | reconstructed fODFs, **healthy**: Lmax down the rows, truth in column 1 then one column per SNR |
| 5 | the same for **edema** |
| 6 | **bias, spread and spurious peak count against SNR** — one curve per Lmax, **one row per kernel**, y limits shared per column |

**Only Figure 1 uses a shared radial scale**, because it is the one panel whose
point is that the edema response is genuinely smaller than the healthy one.
Everywhere else each panel autoscales and fills its box: Figures 2 and 3 are
about the *shape* of the signal at each shell and SNR, Figures 4 and 5 about the
shape of the reconstruction, and shrinking one of them into illegibility would
buy nothing.

Getting that shared scale to work takes two things and only one is obvious.
Scaling a glyph's radius does nothing on its own, because `axis equal` fixes a
panel's aspect ratio but *not* its limits — MATLAB autoscales each subplot to
its own data, so a glyph half the size gets an axis range half as wide and is
drawn at exactly the same size. Figure 1 therefore scales every glyph to a
maximum radius of 1 **and** pins each panel to `GLYPH_LIM`. Its ground-truth
fODF is deliberately *not* on the response scale, since an fODF and a signal are
different quantities.

Set `MAKE_FIGURES = false` in that section to skip them. In MATLAB they appear
inline in the Live Script; under Octave graphics are often unavailable, and the
printed numbers are the deliverable either way — every figure is drawn from the
same arrays Step 7 prints, so the tables and the plots cannot disagree.

**In `smi_manuscript_60deg.m` specifically**, the CSD and MSMT-CSD arms are
still not wired up; what is there is a commented hook just before Step 7. Two
notebooks now run those arms for real — `smi_manuscript_60deg.ipynb` on the
Octave kernel and `csd_manuscript_60deg.ipynb` on the Python kernel — and
**three of the claims that commented hook makes are false**, listed at the top
of this file: the arms do *not* appear automatically in the figures, the peak
finder was *not* scale-free, and the response must be evaluated at MRtrix's own
per-shell b. Take the working notebooks over the stub.

The route into the `.m` file that does work today is
`read_csd_export` + `figures_csd_arms`, which read what the Python notebook
exported. Figure 1 is still the panel built to take an estimated response.

### Why the edema arm looks noisier, when the noise is identical

`sigma = 1/SNR` against `S0 = 1` in both runs, so the noise added is exactly the
same. What shrinks is the **anisotropic signal**. All three edema parameters push
the same way: `f` 0.60 → 0.10 removes most of the stick compartment, which is the
anisotropic one; `Deperp` 0.50 → 1.15 makes the extra-axonal tensor nearly
isotropic; and `fw` 0.02 → 0.35 adds purely isotropic signal. `K_2`, `K_4`, `K_6`
all shrink while `K_0` stays comparable — and the orientation information lives
entirely in those `l >= 2` invariants.

Deconvolution then amplifies rather than rescues, because recovering `p_lm` means
**dividing by `K_l`**: the gain is `1/K_l`, so smaller invariants amplify the same
noise more. Same mechanism as the CSF blow-ups documented elsewhere in this repo
(`g_2 = 1/||K_2||`, 3.6 in CSF against 0.6 in white matter).

Figure 1's shared scale is what makes this visible, and the table it prints gives
the numbers: on that scale the healthy response peaks at 0.83 / 0.74 / 0.68 for
b = 1 / 2 / 3, the edema response at 0.29 / 0.16 / 0.12, and both are 1.000 at
b = 0 by construction. So the framing is not
"edema is noisier" but **"edema has less anisotropic signal, and the
deconvolution amplifies noise in inverse proportion to it."**

## Running it

```
cd deconv_comparison/notebooks
octave-cli --no-gui -q smi_simulation_walkthrough.m      # about 6 minutes
```

or from MATLAB, from anywhere:

```matlab
run('deconv_comparison/notebooks/smi_simulation_walkthrough.m')
```

Nothing needs to be generated first and no data is downloaded. The acquisition
is a **real HCP 3-shell scheme**, read from `../protocol/hcp_real_3shell.txt`
via `mc_config.m`, so the MATLAB side of the package needs no Python.

Almost all the runtime is Step 6. In `smi_simulation_walkthrough.m` that is one
`SMI.fit` per Lmax, and `SMI.fit` is dominated by training its regression rather
than by voxel count, so each of the three Lmax values costs the same ~2 minutes
whether there are 100 voxels or 10,000.

**`smi_manuscript_60deg.m` is a different order of magnitude**, and at its
defaults it runs in hours rather than minutes. It fits
`NKERN × numel(SNR_LIST) × numel(LMAX_LIST)` times — **42** at two kernels,
`SNR_LIST = [5 10 20 30 50 100 Inf]` and `LMAX_LIST = [4 6 8]` — each on `NREP`
voxels, and at
`NREP = 1000` the per-voxel constrained deconvolution is no longer negligible
next to the training. **`NREP` is the knob**: drop it to 25 for a version that
runs while you read it. The file prints the elapsed time of every fit as it
goes, so the first two lines of Step 6 tell you what the whole run will cost.

One entry of the sweep is cheaper than it looks and dearer than it should be:
at `SNR = Inf` there is no noise, so all `NREP` realisations of a condition are
the same signal and come back bit identical. That block costs a full fit to
produce one distinct answer. It earns its place as the reference column of
Figures 3–5 and as the check that the fit reaches the truth at all, but it is
the cheapest entry to drop if runtime is what is hurting.

MRtrix is **not** required. The last step writes the fODFs as MRtrix SH images
so you can check them with `sh2peaks` or `mrview` yourself; the walkthrough
prints those commands and runs none of them.

## Opening it as a MATLAB Live Script

These are plain `.m` files written with publish-style markup (`%%` sections,
`*bold*`, `|monospace|`), so MATLAB converts them with no edits:

- **In the editor**: open the `.m`, then *Save As* and choose `.mlx`.
- **From the command line**:
  ```matlab
  matlab.internal.liveeditor.openAndSave( ...
      'smi_simulation_walkthrough.m', 'smi_simulation_walkthrough.mlx');
  ```

The `.m` is kept as the source of truth rather than the `.mlx` because `.mlx`
is zipped XML: it cannot be diffed, reviewed in a pull request, or run under
Octave, and this package's whole claim is that an outsider can check it.

## What "verifiable" means here

Every step ends with one or more `CHECK` lines that compare its output against
something computed a *different* way, and print `ok` or `** FAILED **` next to
the number. Reading only the `CHECK` lines is a complete audit. The ones worth
understanding:

| check | why it is not circular |
|---|---|
| zonal response == SMI forward model | the same single-fibre signal built from the response-function form and from `SMI`'s own convolution. This identity is what lets an SMI kernel be handed to MRtrix at all |
| harmonics vs direct convolution | the signal computed with no spherical harmonics anywhere, by summing over the quadrature grid. Catches a basis or normalisation error that a self-consistent SH round-trip would hide |
| spherical mean is orientation-blind | SMI fits the kernel from rotational invariants, so the `l = 0` invariant must not depend on where the fibres point |
| recovered sigma | the noise level read back out of the simulated data, rather than the value that was typed in. In the manuscript file this runs once per SNR |
| the ceiling | the same peak finder run on the noise-free truth *truncated to the same Lmax as the fit*. Separates "the method failed" from "this angular order cannot represent the answer" |
| SMI bins the shells | the real protocol has 18 distinct b values; `SMI.Group_dwi_in_shells_b_beta_TE` must recover `[18 90 90 90]`, or every rotational invariant downstream is built from a fraction of the directions |
| the noise-free arm reaches the truth | *(manuscript file only)* with `sigma = 0` the only things left between fit and truth are the band limit and the estimated kernel, so the `SNR = Inf` row must recover the true fibre count wherever the truth itself resolves. If it does not, no finite-SNR row below it is interpretable |
| the noise-free arm is deterministic | *(manuscript files only)* at `sigma = 0` every realisation of a condition is fed the same signal, so the block must come back the same. **See the correction below — the `.m` file states this invariant in a form that is false, and the `.ipynb` states the measured one** |
| more noise does not help | *(manuscript file only)* the mean angular error at the top of the sweep must not exceed the error at the bottom. Compared at the extremes rather than pairwise, so Monte Carlo error between adjacent SNRs cannot trip it |

In the manuscript file the ceiling is not scored by a *copy* of the peak finder:
the band-limited truth is prepended as the first row of every scored block, so
it goes through the same lines as the realisations. That is why the `SNR = Inf`
row and the `ceiling` row of each table print the same angular error.

### Correction: `SMI.fit` is not bit-reproducible voxel to voxel

`smi_manuscript_60deg.m:1069` checks that the noise-free block has standard
deviation **exactly** `0`, and the table above used to describe that as the
invariant. **Both are wrong**, in two separate ways, and the notebook's smoke
run found them. The `.m` file has never been run at full size, so neither had
ever been exercised.

1. **`std(x) == 0` is not a determinism test.** `std` computes `sum(x)/n` first,
   and for *n* identical values that division need not return the value exactly.
   A perfectly deterministic block can yield `std ~ 1e-16`. Whether it does
   depends on the bit pattern of the particular angular error — which is why
   this check failed for two arms and passed for a third on identical logic.

2. **`SMI.fit` genuinely is not bit-reproducible.** Measured directly: feed 27
   voxels on a `[3 3 3]` grid a bit-identical signal at `sigma = 0`, with the
   non-negativity constraint converging in the same 4 iterations in every
   voxel, and SMI returns **three distinct answers — one per slice** — differing
   by `1.6e-12` in `plm` and `1.1e-11` in the kernel. That is BLAS blocking:
   reduction order inside the matrix operations depends on where a voxel sits
   in the array. It is twelve orders of magnitude below the peak finder's 2.6°
   grid resolution and changes no result.

**Both MRtrix arms *are* bit identical**, because `dwi2fod` works one voxel at a
time. So the corrected check is worth keeping rather than deleting: it separates
a method that is exactly reproducible voxel to voxel from one that is
reproducible only to floating point. The notebook tests the SH coefficients
directly against a `1e-9` tolerance and prints the measured number.

Three results in the output are *findings, not failures*, and each is called out
in place because it looks alarming otherwise:

- **The band-limited ground truth is negative over roughly 40% of the sphere.**
  Truncating a Watson mixture at `Lmax` rings, and the rings cross zero. The
  non-negativity constraint is therefore a regularizer, not a statement of
  fact — the truth does not satisfy it either.
- **The 30 degree crossing is resolved 0% of the time, and that is correct.**
  The ceiling shows the noise-free truth also comes back as a single peak: two
  Watson populations at `kappa = 16` merge into one lobe well before any
  deconvolution is involved. This is a property of the *ground truth*, not of
  any method, and it is worth knowing where the boundary sits — measured on the
  truth at Lmax 8, a 30 degree crossing needs `kappa >= 48` (about 8 degrees of
  dispersion) before it separates at all, while 45 and 60 degrees separate at
  every `kappa` from 8 upward. `mc_config.m` keeps `kappa = 16` because real
  white matter disperses; the cost is that the 30 degree condition currently
  measures the band limit rather than the method.
- **A loud `WARNING` about gradient directions on every run.** The supplied
  `.bvec` is unit only to `1.1e-6`, because it was written at seven significant
  figures. `mc_config.m` warns and normalises. It is deliberately noisy: left
  uncorrected it degrades the zonal-response identity at Lmax 8 from `1e-15` to
  `5e-7`, and that failure is invisible unless something checks for it.

## The SNR sweep, and the one trap in it

`smi_manuscript_60deg.m` runs `SNR_LIST = [5 10 20 30 50 Inf]` as six
independent blocks of voxels, one `SMI.fit` per block per Lmax. Each block is
seeded separately (`mc_config.m`'s seed plus the SNR index) so a noise level is
reproducible on its own and adding or removing an entry does not silently
change the realisations of the others. The cost of that choice is that the
blocks do *not* share common random numbers, so a difference between two SNRs
carries the Monte Carlo error of both.

Four numbers are reported per cell, because a sweep is where they stop agreeing:
**correct fibre count** (what tractography consumes), **mean angular error**,
**its standard deviation**, and **spurious peaks per voxel**. A method can hold
its angular error while it starts inventing peaks, or lose the peak entirely and
report a confident wrong direction; no single column catches both.

**"Bias" here is the mean angular error, which is the mean of a non-negative
quantity.** It does not go to zero for a perfect estimator — the floor is set by
the direction grid (about 1.5° at 1500 directions) and by the band limit. Read
it against the noise-free point of the same curve, not against zero. That is
what the `SNR = Inf` column is for.

**The trap: SMI does not fit the kernel at the sigma you pass in.** It
normalises sigma by the measured `b = 0` signal, bins the result into `Nlevels`
equal bins spanning `sigma_norm_limits`, and trains one polynomial regression
per occupied bin, evaluated at the **bin centre** (`SMI.m:2222-2302`). With the
shipped `sigma_norm_limits = [0 0.2]` (`SMI.m:562`) and `Nlevels = 10`
(`SMI.m:388`) that has two consequences only a sweep makes visible:

- **two SNRs can land in the same bin** and then share a regression trained at
  one noise level — σ = 0.02 (SNR 50) and σ = 0.033 (SNR 30) both fall in bin 2;
- **an SNR below 5 is off the end of the trained range.** σ = 0.2 is exactly the
  top edge, so anything noisier is clamped into the top bin and fitted with a
  regression trained at a *lower* noise level than the data actually has.

Step 5 prints the bin each SNR will land in and flags any that are clamped. The
bins are nominal: sigma is normalised by the *measured* `b = 0` signal, so
voxels of one SNR can straddle an edge.

## Why the manuscript uses 60 degrees only

30 and 45 stay in the general walkthrough, but neither makes a clean figure.
Measured on the **noise-free ground truth**, with no fitting and no noise:

| angle | Lmax 4 | Lmax 6 | Lmax 8 | true |
|---|---|---|---|---|
| 30° | 1 | 1 | 1 | 2 |
| 45° | 1 | 2 | 2 | 2 |
| 60° | 2 | 2 | 2 | 2 |

At `kappa = 16` a 30 degree crossing does not separate at any angular order —
it needs `kappa >= 48`, about 8 degrees of dispersion. 45 separates only from
Lmax 6 upward. **60 separates everywhere**, which is what makes it the
configuration where a difference between methods is attributable to the method
rather than to the band limit.

## Relationship to the pipeline

The walkthrough is not a copy of `gen_montecarlo.m`; it reads its constants
from `../mc_config.m`, the same file the pipeline reads, so the two cannot
drift apart about what is being simulated. The SNR sweep keeps that
correspondence: `gen_montecarlo.m` takes one SNR per invocation, and the
manuscript file's Step 6 makes one `SMI.fit` call per SNR for the same reason —
no fit ever sees two noise levels at once. It does differ from the pipeline in
two deliberate ways: peaks are found on a fixed direction grid so the operation
is visible in the file, where the pipeline uses `sh2peaks` for every arm; and it
fits at three angular orders where the pipeline fixes Lmax 6 to keep the three
methods comparable. Step 8 exports the fODFs so `sh2peaks` can be run on them
directly.

The manuscript file's export covers the **whole sweep** in one image per Lmax —
every SNR block, in the order Step 5 laid them out — so a single `sh2peaks` run
reaches all noise levels. `export/voxel_key.txt` therefore carries an **SNR
column** alongside the crossing angle and the true fibre axes, and its header
records how many contiguous blocks of how many voxels are in the file, and in
what order. The noise-free arm appears there as `Inf`.

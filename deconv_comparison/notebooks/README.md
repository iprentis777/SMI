# Notebooks

Narrated, step-by-step versions of the simulations, for reading rather than
running in bulk. They are the answer to "how do I check this myself?".

| notebook | what it walks through |
|---|---|
| `smi_simulation_walkthrough.m` | the SMI arm end to end on a real HCP protocol: ground truth, kernel, forward convolution, Rician noise, `SMI.fit` at Lmax 4/6/8, peaks, MRtrix export. Conditions: single fibre, 30, 45, 60 degrees. **One SNR** |
| `smi_manuscript_60deg.m` | the manuscript cut — **single fibre and 60 degrees only**, **swept across SNR** — with six isometric figures |

`smi_simulation_walkthrough.m` draws four figures: the ground truth fibre
geometry, the kernel as a response function, the forward convolution with and
without noise, and the reconstructed fODFs at each Lmax.

`smi_manuscript_60deg.m` draws six, and every 3D panel opens in an **isometric**
view so all subfigures are readable without touching the camera:

| figure | layout |
|---|---|
| 1 | the two ground truth fibre configurations |
| 2 | a montage of response glyphs, b = 0 to 3 |
| 3 | the signal as surfaces on the sphere: **b down the rows**, SNR across the columns, increasing left to right and ending at the ground truth |
| 4 | reconstructed fODFs, **single fibre**: Lmax down the rows, the ground truth in the first column and then one column per SNR |
| 5 | reconstructed fODFs, **60 degree crossing**, same layout |
| 6 | **bias, spread and spurious peak count against SNR** — one curve per Lmax, one row per condition |

Set `MAKE_FIGURES = false` in that section to skip them. In MATLAB they appear
inline in the Live Script; under Octave graphics are often unavailable, and the
printed numbers are the deliverable either way — every figure is drawn from the
same arrays Step 7 prints, so the tables and the plots cannot disagree.

The CSD and MSMT-CSD arms are not covered yet. Figure 2 is the panel built to
take their responses and Figures 4–6 the panels built to take their fODFs.

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
`numel(SNR_LIST) × numel(LMAX_LIST)` times — 18 at `SNR_LIST = [5 10 20 30 50
Inf]` and `LMAX_LIST = [4 6 8]` — each on `2 × NREP` voxels, and at
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
| the noise-free arm is deterministic | *(manuscript file only)* at `sigma = 0` every realisation is the same signal and `SMI.fit` trains one regression per call, so the whole block must come back bit identical — its standard deviation is checked to be exactly `0` |
| more noise does not help | *(manuscript file only)* the mean angular error at the top of the sweep must not exceed the error at the bottom. Compared at the extremes rather than pairwise, so Monte Carlo error between adjacent SNRs cannot trip it |

In the manuscript file the ceiling is not scored by a *copy* of the peak finder:
the band-limited truth is prepended as the first row of every scored block, so
it goes through the same lines as the realisations. That is why the `SNR = Inf`
row and the `ceiling` row of each table print the same angular error.

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

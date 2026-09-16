# The simulation

One file: [`smi_wm_60deg.m`](smi_wm_60deg.m).

It is written to be read while it runs. Plain `.m`, so it runs in MATLAB and GNU
Octave and converts to a Live Script unedited. **Its Step 0 documents every
design decision and the measurement behind it** — read that rather than
duplicating it here.

## What it asks

> Given the same response function and the same data, how do the deconvolution
> algorithms compare?

Everything that made that question unanswerable in earlier versions was removed
rather than argued about. See
[`Archive/old_simulations/README.md`](../../Archive/old_simulations/README.md)
for the three defects that forced the rewrite.

## The arms

| arm | what it is |
|---|---|
| **SMI fixed** | `SMI.get_plm_from_S_and_kernel` with a **fixed** kernel — the same one the CSD arms get as a response. No Rician bias correction, so the same noise model as MRtrix |
| **SMI fitted** | `SMI.fit`, estimating the kernel per voxel. Off by default (`RUN_ARM2 = 0`) because it is expensive: one `SMI.fit` per Lmax per SNR |
| **SSST-CSD** | `dwi2fod csd` on the top shell |
| **MSMT def** | `dwi2fod msmt_csd` at MRtrix's shipped `-neg_lambda 1e-10 -norm_lambda 1e-10` |
| **MSMT tuned** | the same, at `-neg_lambda 1 -norm_lambda 1e-3` |

The difference between the two SMI arms is **the cost of estimating the kernel**,
measured against a fair baseline rather than asserted.

Both MSMT variants run because the difference between them is **order
dependent** — at Lmax 4 the defaults fail outright (0.0% of crossings resolved),
at Lmax 6 they are the most noise-robust arm (83.9% correct at SNR 5 against
SMI's 44.1%) at the cost of high-SNR bias. Report MSMT numbers with both lambda
values *and* the Lmax; none of the three generalises to the others.

## Reading it

**Every step ends in one or more `CHECK` lines** comparing its output against
something computed a different way. Reading only the `CHECK` lines is a complete
audit. They are the pass/fail signal — two worth watching:

- Step 6a, *"unconstrained fixed-kernel fit inverts the forward model"* — must be
  `< 1e-12`. If this fails the kernel-to-response conversion or the design matrix
  is wrong and no number in the file is trustworthy.
- Step 7, *"the ceiling is orientation invariant"* — the band-limit ceiling must
  not depend on where the crossing points. This is the check that would have
  caught the old fixed-grid peak finder.

## Running it

```sh
cd deconv_comparison/notebooks
octave-cli --no-gui -q smi_wm_60deg.m
```

⚠️ **It defaults to the full sweep and runs for hours** — 18 orientations × 50
reps × 16 SNR at Lmax 4, 6 and 8. Set `SMOKE_TEST = true` for a reduced run in
minutes with every `CHECK` still executed and indicative numbers only.

The statistical argument for the full size is in the file's own `SNR_LIST`
comment: 900 voxels per SNR give `correct(%)` a worst-case standard error of
**1.7 percentage points**, against **7.2** at smoke-test size. Steps in a curve
smaller than ~3.4 pp are sampling noise, not structure.

**Requires MRtrix3 on the PATH** (`dwi2fod`, `sh2peaks`, `mrinfo`, `dwiextract`)
or it stops at Step 6c. Nothing here reimplements MRtrix; the only MRtrix
behaviour implemented locally is reading and writing its image format
(`../mrtrix_io.m`), which `mrinfo` checks on every run.

## Conventions that will silently ruin a result

- **`CS_phase = 0`.** At SMI's default of 1 the SH basis differs from MRtrix's by
  `(-1)^m` — a 180° rotation about z of every fODF. A basis error shows as ~71°
  on every arm at once.
- **`p_00 = 1`.** SMI's fODF is normalized; an MRtrix FOD is not, and its `l = 0`
  varies per voxel. The peak threshold subtracts each voxel's *own* `l = 0`, not
  the constant `1/(4π)` — subtracting the constant would silently mis-threshold
  every CSD arm.
- **The band-limited truth is itself negative** over part of the sphere.
  Non-negativity is a regularizer, not a statement of fact; the truth does not
  satisfy it either.

## Historical documentation

The 699-line notebook README that covered the earlier simulations — MRtrix
convention findings, the CSD response derivation, and a correction on `SMI.fit`
not being bit-reproducible voxel to voxel — is preserved at
[`Archive/old_simulations/NOTEBOOKS_README.md`](../../Archive/old_simulations/NOTEBOOKS_README.md).

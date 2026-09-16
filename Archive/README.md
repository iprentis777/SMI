# Archive

This directory indexes exploratory work and superseded workflows that are not
part of the recommended analysis pipeline. Material is retained when it remains
useful for understanding SMI, reproducing an older result, or explaining why a
direction was not pursued.

Archive status is not a claim that the work is incorrect. It means that the
exercise is primarily useful for understanding the model, has not been
validated enough for routine use, or no longer reflects the direction of the
main pipeline.

## Superseded workflow: original deconvolution comparison

[`deconv_pipeline/`](deconv_pipeline/) contains the original Octave, MRtrix,
and Python Monte Carlo campaign. It produced
[`Reports/deconv_tables.md`](../Reports/deconv_tables.md) and
[`Reports/REPORT_SMI_deconvolution_MonteCarlo.md`](../Reports/REPORT_SMI_deconvolution_MonteCarlo.md),
so it remains the provenance for those numbers.

It is not the active comparison and its paths need repair before it can run from
its archived location. New work should use
[`deconv_comparison/notebooks/smi_wm_60deg.m`](../deconv_comparison/notebooks/smi_wm_60deg.m),
which runs every arm on the same simulated data in one script. See the
[`deconv_pipeline` README](deconv_pipeline/README.md) for the exact boundary.

## Retired: anisotropy modulation of the fODF

[`patch_history/fODF_modulation/`](patch_history/fODF_modulation/) holds the
removed implementation, its measurement report and its example.

SMI stores a normalized fODF with `p_00 = 1`, so its isotropic floor is a fixed
`1/(4*pi) = 0.0796` — **above** MRtrix's default iFOD2 `-cutoff` of 0.05. An
unmodulated fODF therefore passes the tractography termination test everywhere
in the brain, CSF included. Anisotropy modulation multiplied the fODF by a
per-voxel coherence weight to restore that amplitude information without keying
on tissue type, which would delete edema.

This was previously *semi-retired* — off by default but still shipped. It is now
removed: `options.fODF_modulation` is no longer read by any code path.

It is retired because its evidence is simulation-only and was never validated on
real edema, its preferred `p2product` weight fails for some symmetric fibre
geometries, it changes the coefficient convention in density mode, and it cannot
stabilize an ill-conditioned deconvolution — the weight is clipped at 1, so
`w * 1e13` is still `1e13`. Regularization is what prevents that.

The exercise established several useful negative results, which are the reason
it is archived rather than deleted: `p4` is not a reliable modulation weight and
is *worse than no weighting at all*; tissue-fraction weights suppress the edema
class they were intended to preserve; the original high-order loss was not
caused by Tikhonov damping but by the non-negativity constraint; and
`degenerate = 'clip'` was a bad default that gave blown-up voxels the maximum
weight in the volume.

Two things deliberately stayed behind. `SMI.grab_pl` and `SMI.grab_kernel_pl`
are kept — modulation was their only caller, but they are generic public
accessors. And `helpers/fODF_modulation_helpers.m` was **renamed** to
[`helpers/fODF_sim_helpers.m`](../helpers/fODF_sim_helpers.m) rather than
archived: despite the old name it contains no modulation code, only the generic
forward-simulation toolkit that eleven active files depend on.

## Retired: Tikhonov damping of the fODF deconvolution

[`patch_history/fODF_tikhonov/`](patch_history/fODF_tikhonov/) holds the removed
implementation and the measurements that removed it.

The fODF deconvolution is ill-conditioned, and Tikhonov damping added a penalty
`lambda_tikhonov^2*||Gamma*plm||^2` intended to suppress the high-order
coefficients the kernel attenuates most. It was off by default and is now gone
entirely.

It is retired because two independent measurements found it inert: swept from 0
to 0.8 at three noise levels the high-`l` bands were identical to three
decimals, and `0.3` against `0` moved a 45 degree error from 21.28 to 21.27
degrees. The first of those was made while investigating why `l = 4` power was
being lost, and it overturned the explanation then in the code — the cause is
the non-negativity constraint together with error in the estimated kernel, not
damping.

The useful negative result is that **non-negativity is the only fODF regularizer
in SMI that demonstrably changes a result**, so it is the only one worth tuning.
Note also that removing this leaves SMI less regularized than the CSD arms it is
compared against, since `dwi2fod csd` ships `-norm_lambda 1` — Tikhonov with
`Gamma = I`.

A caller that still sets `lambda_tikhonov` is **silently ignored**, not
rejected. See the directory README before re-running an old analysis against
old numbers.

## Learning exercise: viewing the response kernel as zonal harmonics

This exercise expresses the fitted SMI kernel as the zonal harmonic response
used by CSD tools. For a single fibre along `z`,

```text
R(theta) = sum_l K_l(b) (2l+1) P_l(cos theta)
         = sum_l r_l Y_l0(theta)

r_l = K_l(b) * sqrt((2l+1) * 4*pi)
```

It was valuable for establishing conventions, checking the SMI forward model,
and understanding how a parametric SMI kernel relates to an MRtrix response
file. It is archived because it is explanatory material rather than a step a
reader needs in order to run the toolbox.

Artifacts:

- `examples/example_SMI_response_shview.m` — profiles, glyphs, compartment
  decomposition, and MRtrix response export
- `helpers/SMI_response_helpers.m` — kernel/response conversion helpers
- `tests/test_SMI_response_helpers.m` — eight convention and round-trip checks
- `Figures/fodf_response_shview.png` — generated figure

The example accepts a kernel from a real fit through
`SMI_response_helpers().kernel_from_out`. It also checks its zonal
reconstruction against SMI's own forward model before drawing anything.

## Superseded simulations

[`old_simulations/`](old_simulations/) holds every simulation built before
`smi_wm_60deg.m`, including the previous manuscript source and the step-by-step
walkthrough. Its README records the three measured defects that forced the
rewrite: a fixed-grid peak finder that added +/-1.5 deg of orientation-dependent
noise and inverted an arm ranking, a single crossing orientation that nearly
doubled an apparent gap between methods, and arms that were recovering
*different objects* because only one of them got a dispersion-matched response.

## Patch history

[`patch_history/`](patch_history/) holds the raw development history as patches,
plus one subfolder per removed feature. New contributors do not need it; it
exists so every number in the reports has a provenance and so a removed feature
can be recovered rather than only described.

## Adding future exercises

When an exploratory feature leaves the active pipeline, add it here with:

1. its status (`learning exercise`, `semi-retired`, or `retired`);
2. the question it was meant to answer;
3. what was learned, including negative results;
4. why it is not recommended pipeline guidance; and
5. links to the implementation, tests, reports, and figures that preserve its
   provenance.

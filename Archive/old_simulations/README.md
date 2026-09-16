# Superseded simulations

Every simulation this fork built before `smi_wm_60deg.m`, kept because the
reports cite them and because several carry findings with no other provenance.

**None of these is the active simulation.** That is
[`deconv_comparison/notebooks/smi_wm_60deg.m`](../../deconv_comparison/notebooks/smi_wm_60deg.m),
which asks one question — *given the same response function and the same data,
how do the deconvolution algorithms compare?* — and removes everything that made
that question unanswerable.

## What is here

| file | what it was | why it is here |
|---|---|---|
| `smi_manuscript_60deg.m` | the previous manuscript source. 60° crossing, two kernels (healthy + edema), SNR swept, all three arms in one file | **directly superseded by `smi_wm_60deg.m`.** See below for the three defects that forced the rewrite |
| `smi_simulation_walkthrough.m` | the SMI arm taken apart step by step on a real HCP protocol. Single fibre plus 30/45/60°, one SNR, Lmax 4/6/8 | the teaching version. Its `CHECK`-per-step structure is what `smi_wm_60deg.m` inherited |
| `sweep_deconv_settings.m` | one-at-a-time and grid sweeps over the deconvolution settings | superseded as a tuning tool; `lambda_nonneg` is the only weight left to tune |
| `smi_free_l0_experiment.m` | the experimental free `l = 0` deconvolution, where `p_00` is estimated rather than fixed at 1 | pairs with `SMI_freeL0.m` and `helpers/fODF_free_l0_deconv.m`, **both of which are still in the active tree** |
| `check_manuscript_static.m` | static checks on `smi_manuscript_60deg.m` — that it parses, that scoring arrays are subscripted correctly | bound to the file it checks |
| `measure_glyph_spread.m` | how much drawn glyph radius varies between noise realisations vs across SNR | measures a figure of the archived file |
| `NOTEBOOKS_README.md` | the 699-line notebook documentation | **read this before reusing anything here.** It holds the MRtrix convention findings, the CSD response derivation, and the non-reproducibility correction |

## Why `smi_manuscript_60deg.m` was replaced

Three defects, each measured rather than argued:

1. **Peaks came from a fixed 1500-direction grid with no refinement.** On the
   noise-free band-limited truth at Lmax 6 that scorer reports the ceiling as
   anything from **0.488 to 2.847 deg** depending on where the crossing sits
   relative to the grid, where `sh2peaks` reports **1.300 deg at every
   orientation** (spread 0.010 deg over 18). The grid was adding ±1.5 deg of
   orientation-dependent noise to every angular number, and it inverted at least
   one arm ranking.
2. **The crossing was simulated at one orientation.** At SNR 10 that single
   orientation was simultaneously near-best for SMI and exactly worst for
   SSST-CSD, nearly doubling the apparent gap between them.
3. **The arms were not recovering the same object.** CSD got a
   dispersion-matched response while SMI deconvolved with the delta kernel, so
   the shared "ceiling" was wrong for the CSD arms at Lmax 8 — 1.27 deg printed
   against their true 2.82.

`smi_wm_60deg.m` fixes all three: `sh2peaks` for every arm including SMI's, 18
orientations, and one shared response handed to every arm.

## Running them

They will not run unmodified. Beyond the usual path repair, two toolbox features
they use **no longer exist**:

- `options.fODF_regularization.lambda_tikhonov` — removed. It is *silently
  ignored*, not rejected, so these files will run and quietly differ. See
  [`../patch_history/fODF_tikhonov/`](../patch_history/fODF_tikhonov/).
- `options.fODF_modulation` and `SMI.modulate_fODF` — removed. See
  [`../patch_history/fODF_modulation/`](../patch_history/fODF_modulation/).

Their calls to the simulation helpers were updated when
`helpers/fODF_modulation_helpers.m` was renamed to `helpers/fODF_sim_helpers.m`,
so that part is current.

## What stayed active for these files' sake

`SMI_freeL0.m` and `helpers/fODF_free_l0_deconv.m` are still in the active tree
even though their only in-repo driver (`smi_free_l0_experiment.m`) is archived
here. They were left alone deliberately rather than archived by association —
the free `l = 0` question is independent of this cleanup.

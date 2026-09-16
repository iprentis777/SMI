# Additions to the NYU SMI toolbox

This fork adds two opt-in features to the fODF deconvolution. Both are **off by
default**, so a fit that does not ask for them is bit-identical to stock NYU SMI
(verified, difference exactly 0).

| addition | what it does | default | where |
|---|---|---|---|
| [Regularized deconvolution](fODF_regularization/) | a non-negativity constraint on the fODF, as in constrained spherical deconvolution | off (`flag_nonneg = 0`) | `options.fODF_regularization` |
| [Outlier capping](fODF_outlier_cap/) | reduces isolated pathologically bright glyphs, without touching orientation | off (`flag_cap = 0`) | `options.fODF_outlier` |

Everything else in this repository is either the original toolbox, the
simulation that measures these two features
([`deconv_comparison/`](../deconv_comparison/)), or archived material
([`Archive/`](../Archive/)).

## The shortest possible version

```matlab
options.flag_fit_fODF = 1;

% Addition 1: non-negativity constraint
options.fODF_regularization.flag_nonneg = 1;

% Addition 2: post hoc outlier cap
options.fODF_outlier.flag_cap = 1;

out = SMI.fit(dwi, options);
```

Both write their settings and per-voxel diagnostics back into `out`, and both
leave `out.plm` in the original convention — the corrected fODF comes back
separately as `out.plm_capped`. See each subfolder for what the numbers mean and
what was measured.

## Why these two and not others

Two further additions were built, measured, and **removed**. They are documented
in [`Archive/patch_history/`](../Archive/patch_history/) rather than deleted,
because the measurements that retired them are useful:

- **Tikhonov damping** of the deconvolution — inert at the weights it shipped
  with, and actively harmful above them.
- **Anisotropy modulation** of the fODF — addressed a real problem (a normalized
  fODF's isotropic floor sits above MRtrix's tractography cutoff) but was never
  validated outside simulation.

The pattern is deliberate: an addition stays only while a measurement says it
earns its place.

## What is coming

This folder will also carry **curated patch files** showing exactly what changes
between the NYU implementation and this one, feature by feature. Until then, the
raw development history — every patch, in order — is in
[`Archive/patch_history/`](../Archive/patch_history/), and the root
[`README.md`](../README.md) documents both features in full.

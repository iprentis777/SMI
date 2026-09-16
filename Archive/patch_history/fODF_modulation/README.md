# Retired: anisotropy modulation of the fODF

**Status:** retired. Removed from the toolbox; `options.fODF_modulation` is no
longer read by any code path.

Previously this was *semi-retired* — off by default but still shipped. It is now
archived outright.

## The question it was meant to answer

SMI stores a normalized fODF with `p_00 = 1`, so every voxel has the same total
fODF mass. That has a concrete downstream consequence:

```text
isotropic floor of a normalized fODF = 1/(4*pi) = 0.0796
MRtrix iFOD2 default -cutoff        = 0.05
```

The floor is **above** the cutoff, so an unmodulated SMI fODF passes the
tractography termination test *everywhere in the brain, CSF included* — 100% of
simulated GM and CSF voxels survived it. Modulation multiplied the fODF by a
per-voxel orientational-coherence weight so that amplitude could separate
coherent white matter from isotropic tissue, restoring the information the
normalization removes.

Weighting by *coherence* rather than tissue type was deliberate, and was the
idea's strongest feature: fibres displaced by edema stay coherent and keep a
high weight even when their axonal fraction has collapsed, whereas an `f` or
`1-fw` threshold deletes them. Measured, `f` kept 94% of a simulated edema class
and `1-fw` kept 85%, against **100%** for every `p2` weight.

## Why it is retired

Not because it was measured wrong — the negative results below are real and
worth keeping. It is retired because the repository no longer needs it, and
because it carried costs that were never paid off:

- its evidence is **simulation-only** and was never validated on real edema;
- its preferred `p2product` weight **fails for some symmetric fibre geometries**;
- in `density` mode it **changes the coefficient convention** by including and
  rescaling the `l = 0` term, so a modulated fODF is not in the convention the
  rest of the toolbox assumes;
- its output basis and downstream MRtrix use require explicit convention checks;
- modulation **cannot stabilize an ill-conditioned deconvolution** and does not
  replace regularization — the weight is clipped at 1, so `w * 1e13` is still
  `1e13`;
- a reviewer would reasonably ask why this weighting belongs in the pipeline
  instead of a more conventional CSD-style treatment.

## What was learned — the negative results are the valuable part

1. **`p4` is not a usable modulation weight**, in either flavour. The kernel
   `p4` is biased upward, and the deconvolved `pl4` comes back at roughly a
   third of its true value (median 0.26 in single-fibre WM against a true 0.71).
   Weighting by it measures the regularizer rather than the tissue and leaves
   **84% of CSF above cutoff — worse than not weighting at all.**
2. **The `pl4` loss is a non-negativity effect, not a damping effect.** An
   earlier version of this work blamed Tikhonov damping. That was measured wrong
   and is one of the three findings that later retired damping entirely
   (`../fODF_tikhonov/`). The cause is the non-negativity constraint together
   with error in the estimated kernel, whose `K_l` at high `l` is small and very
   sensitive.
3. **Tissue-fraction weights suppress the class they were meant to preserve.**
   See the edema numbers above.
4. **The two `p2` estimates fail in opposite ways and their product cancels
   both.** `kernel_p2` is stable but does not fall to 0 in CSF (median 0.31
   against a true 0), because its polynomial regression returns roughly the
   prior mean where the `l = 2` signal is uninformative. `pl2` is unbiased in
   CSF (median 0.076) but inherits the deconvolution's noise and has a heavy
   tail. Alone they leave 54–80% and 20–41% of CSF above cutoff respectively;
   the product was the only weight that worked.
5. **`degenerate = 'clip'` was a bad default.** In a blown-up voxel the raw `p`
   exceeds the clip, so the voxel receives weight exactly 1.0 — the *maximum* in
   the volume. `'reject'` existed and should have been the default.

## What stays behind

`SMI.grab_pl` and `SMI.grab_kernel_pl` are **kept** in the toolbox. Modulation
was their only in-toolbox caller, but they are generic public accessors and
`grab_kernel_pl` encodes the kernel's TE-dependent index layout, which is useful
independently.

`helpers/fODF_modulation_helpers.m` has been **renamed to
`helpers/fODF_sim_helpers.m`**. Despite its old name it contains no modulation
code at all — it is the generic forward-simulation toolkit (`dirs`, `watson`,
`watson_plm`, `mixture_plm`, `signal`, `peak`) used by eleven files including
`smi_wm_60deg.m` and the test suite. Archiving it by name would have broken the
active simulation.

## Provenance

| file | what |
|---|---|
| `0001-remove-fodf-modulation-solver.patch` | the removal from `SMI.m` and `SMI_freeL0.m`: `fODF_ModulationDefaults`, `fODF_ModulationWeight`, `modulate_fODF`, and all `SMI.fit` wiring |
| `REPORT_fODF_modulation.md` | the full measurement report, moved here from `Reports/` |
| `example_fODF_modulation.m` | the seven-class simulation, moved here from `examples/` |

To restore the implementation, apply the patch in reverse from the repository
root:

```sh
git apply -R Archive/patch_history/fODF_modulation/0001-remove-fodf-modulation-solver.patch
```

**`example_fODF_modulation.m` does not run as archived.** It calls
`SMI.fODF_ModulationDefaults` and `SMI.modulate_fODF`, which no longer exist, so
it needs the patch reverse-applied first. Its helper calls were updated to the
renamed `fODF_sim_helpers`, so that part is current.

## Before restoring it

Revisit this only if a concrete real-data or manuscript question requires it.
Any reactivation should begin with a known edema ROI rather than simulation,
explicit MRtrix basis validation, `degenerate = 'reject'` rather than `'clip'`,
and a reviewer-facing justification for departing from established
regularization and tractography conventions. The underlying problem it
addressed — that a `p_00 = 1` fODF has an isotropic floor above MRtrix's default
cutoff — is real and still unsolved; modulation was one answer to it, not the
only possible one.

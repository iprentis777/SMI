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
[`deconv_comparison/notebooks/smi_manuscript_60deg.m`](../deconv_comparison/notebooks/smi_manuscript_60deg.m),
which runs all three arms on the same simulated data in one script. See the
[`deconv_pipeline` README](deconv_pipeline/README.md) for the exact boundary.

## Semi-retired: anisotropy modulation of the fODF

SMI stores a normalized fODF with `p_00 = 1`, so every voxel has the same total
fODF mass. Anisotropy modulation explored multiplying that fODF by a per-voxel
coherence weight so amplitude could help distinguish coherent white matter from
isotropic tissue during tractography.

The implementation remains available through `options.fODF_modulation` and
`SMI.modulate_fODF`, but it is off by default and is not part of the recommended
pipeline. Treat it as semi-retired because:

- its evidence is simulation-only and has not been validated in real edema;
- its preferred `p2product` weight fails for some symmetric fibre geometries;
- it changes the coefficient convention in density mode by including and
  rescaling the `l=0` term;
- its output basis and downstream MRtrix use require explicit convention checks;
- modulation cannot stabilize an ill-conditioned deconvolution and does not
  replace regularization;
- a reviewer would reasonably ask why this additional weighting belongs in the
  pipeline instead of a more conventional CSD-style treatment.

The exercise did establish several useful negative results: `p4` is not a
reliable modulation weight, tissue-fraction weights can suppress the edema
class they were intended to preserve, and the original high-order loss was not
caused by Tikhonov damping.

Artifacts:

- `examples/example_fODF_modulation.m` — seven-class simulation
- `helpers/fODF_modulation_helpers.m` — reusable simulation helpers
- `Reports/REPORT_fODF_modulation.md` — measurements and limitations
- `SMI.m` — retained opt-in implementation

Revisit this work only if a concrete real-data or manuscript question requires
it. Any reactivation should begin with a known edema ROI, explicit MRtrix basis
validation, and a reviewer-facing justification for departing from established
regularization and tractography conventions.

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

## Adding future exercises

When an exploratory feature leaves the active pipeline, add it here with:

1. its status (`learning exercise`, `semi-retired`, or `retired`);
2. the question it was meant to answer;
3. what was learned, including negative results;
4. why it is not recommended pipeline guidance; and
5. links to the implementation, tests, reports, and figures that preserve its
   provenance.

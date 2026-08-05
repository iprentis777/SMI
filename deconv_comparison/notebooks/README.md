# Notebooks

Narrated, step-by-step versions of the simulations, for reading rather than
running in bulk. They are the answer to "how do I check this myself?".

| notebook | what it walks through |
|---|---|
| `smi_simulation_walkthrough.m` | the SMI arm end to end on a real HCP protocol: ground truth, kernel, forward convolution, Rician noise, `SMI.fit` at Lmax 4/6/8, peaks, MRtrix export. Conditions: single fibre, 30, 45, 60 degrees |
| `smi_manuscript_60deg.m` | the manuscript cut of the same run — **single fibre and 60 degrees only** — with four isometric figures |

Both draw four figures — the ground truth fibre geometry, the kernel as a
response function, the forward convolution with and without noise, and the
reconstructed fODFs at each Lmax. In `smi_manuscript_60deg.m` every 3D panel
opens in an **isometric** view so all subfigures are readable without touching
the camera, Figure 2 is a montage of response glyphs from b = 0 to 3, and
Figure 3 renders the signal as surfaces on the sphere rather than as scatter
plots. Set `MAKE_FIGURES = false` in that section to
skip them. In MATLAB they appear inline in the Live Script; under Octave
graphics are often unavailable, and the printed numbers are the deliverable
either way.

The CSD and MSMT-CSD arms are not covered yet.

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
via `mc_config.m`, so the MATLAB side of the package needs no Python. Almost
all the runtime is Step 6: `SMI.fit` is dominated by training its regression
rather than by voxel count, so each of the three Lmax values costs the same
~2 minutes whether there are 100 voxels or 10,000.

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
| recovered sigma | the noise level read back out of the simulated data, rather than the value that was typed in |
| the ceiling | the same peak finder run on the noise-free truth *truncated to the same Lmax as the fit*. Separates "the method failed" from "this angular order cannot represent the answer" |
| SMI bins the shells | the real protocol has 18 distinct b values; `SMI.Group_dwi_in_shells_b_beta_TE` must recover `[18 90 90 90]`, or every rotational invariant downstream is built from a fraction of the directions |

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
drift apart about what is being simulated. It does differ from the pipeline in
two deliberate ways: peaks are found on a fixed direction grid so the operation
is visible in the file, where the pipeline uses `sh2peaks` for every arm; and it
fits at three angular orders where the pipeline fixes Lmax 6 to keep the three
methods comparable. Step 8 exports the fODFs so `sh2peaks` can be run on them
directly.

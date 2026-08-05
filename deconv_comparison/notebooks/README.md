# Notebooks

Narrated, step-by-step versions of the simulations, for reading rather than
running in bulk. They are the answer to "how do I check this myself?".

| notebook | what it walks through |
|---|---|
| `smi_simulation_walkthrough.m` | the SMI arm end to end: protocol, ground truth, kernel, forward convolution, Rician noise, `SMI.fit`, peaks |

The CSD and MSMT-CSD arms are not covered yet.

## Running it

```
cd deconv_comparison/notebooks
octave-cli --no-gui -q smi_simulation_walkthrough.m      # about 2 minutes
```

or from MATLAB, from anywhere:

```matlab
run('deconv_comparison/notebooks/smi_simulation_walkthrough.m')
```

Nothing needs to be generated first and no data is downloaded. The gradient
table is read from `../protocol/hcp_like_3shell.txt`, which is tracked, so the
MATLAB side of the package needs no Python. The last step calls MRtrix's
`sh2peaks` if it is installed and skips itself if it is not.

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
| the ceiling | the same peak finder run on the noise-free band-limited truth. Separates "the method failed" from "nothing could have succeeded" |
| `sh2peaks` cross-check | the walkthrough's own peak finder against the one every number in the report actually comes from |

Two results in the output are *findings, not failures*, and both are called out
in place because they look alarming otherwise:

- **The band-limited ground truth is negative over roughly 40% of the sphere.**
  Truncating a Watson mixture at `Lmax` rings, and the rings cross zero. The
  non-negativity constraint is therefore a regularizer, not a statement of
  fact — the truth does not satisfy it either.
- **The 15 degree crossing is resolved 0% of the time, and that is correct.**
  The ceiling shows the noise-free truth also comes back as a single peak. Two
  Watson populations with `kappa = 16` merge well before any deconvolution is
  involved.

## Relationship to the pipeline

The walkthrough is not a copy of `gen_montecarlo.m`; it reads its constants
from `../mc_config.m`, the same file the pipeline reads, so the two cannot
drift apart about what is being simulated. It does differ from the pipeline in
one deliberate way: peaks are found on a fixed direction grid so the operation
is visible in the file, where the pipeline uses `sh2peaks` for every arm. Step
8 measures the gap between the two.

# Constrained SMI vs SSST-CSD vs MSMT-CSD

The simulation behind the manuscript, in the design of Jeurissen et al. (2014):
synthesise noise-free signal vectors for crossing white matter fibres by forward
convolution, add complex Gaussian noise so the magnitude is Rician, and
deconvolve many independent realisations of each condition with each method.

**Everything here serves one file**, `notebooks/smi_manuscript_60deg.m`. It runs
all three arms — SMI, SSST-CSD and MSMT-CSD — on **one simulation, one noise
draw, one peak finder, in one scope**. There is no second forward model to keep
in step and no question about whether the arms saw the same data: they are the
same array.

**CSD and MSMT-CSD are MRtrix3 3.0.4 itself.** `dwi2fod`, `dwiextract` and
`mrinfo` are the binaries, called as subprocesses. The only MRtrix behaviour
implemented locally is reading and writing its image format (`mrtrix_io.m`),
which `mrinfo` checks on every run.

**Everything here is simulation.** No result in this directory has touched a
patient scan.

## What is compared

| arm | what it is |
|---|---|
| SMI constrained | `SMI.fit`, `fODF_regularization.flag_nonneg = 1`, `CS_phase = 0`. **Estimates its kernel per voxel** |
| SSST-CSD | `dwi2fod csd` on the top shell |
| MSMT-CSD | `dwi2fod msmt_csd` on all four shells, 3 tissues, `-neg_lambda 1` |

Two settings decide what "fair" means here, and both are documented at length in
`Reports/REPORT_CSD_response_derivation.md`:

- **`CSD_RESPONSE_KERNEL = 'healthy'`** — both tissues are deconvolved with the
  healthy WM response, because that is what a population-averaged single-fibre
  estimate is. This is the one place the arms are not on equal footing, and the
  inequality is real: SMI re-estimates its kernel per voxel and adapts, CSD
  cannot.
- **`MSMT_NEG_LAMBDA = 1`** — MRtrix ships `1e-10`, which leaves `msmt_csd`
  effectively unconstrained while `csd`'s constraint is at strength 1. Running
  both "at their defaults" is not a like-for-like comparison and produced a
  wrong result once already.

**`CS_phase = 0` matters.** At SMI's default of 1 the SH basis differs from
MRtrix's by `(-1)^m`, a 180 degree rotation about z of every fODF.

## The protocol

**A real HCP 3-shell acquisition** — 288 volumes, 18 at b = 5 s/mm² plus 90
directions each at nominal b = 1, 2, 3 ms/µm², read from
`protocol/hcp_real_3shell.txt` through `mc_config.m` so no arm can disagree
about what was acquired.

Three properties of the real scheme, all checked rather than assumed:

- the b = 0 volumes are **b = 5 s/mm², not 0**, and carry unit direction vectors
  even though the direction is meaningless there. `B0_SNAP` sets them to
  exactly 0 so `S(0)/S0 = 1` holds exactly;
- the b values **jitter within each shell** (18 distinct values), and
  `SMI.Group_dwi_in_shells_b_beta_TE` bins them — checked to recover
  `[18 90 90 90]`;
- the supplied `.bvec` is unit only to **1.1e-6**. `mc_config.m` warns loudly
  and normalises. Left alone this breaks any calculation reading `g(3)` as
  `cos θ`: at Lmax 8 the zonal-response identity degrades from 1e-15 to 5e-7.

## What is here

| file | what |
|---|---|
| `notebooks/smi_manuscript_60deg.m` | **the manuscript simulation.** All three arms, two kernels, SNR swept, seven figures |
| `notebooks/smi_simulation_walkthrough.m` | the SMI arm taken apart step by step, one SNR, 30/45/60 degrees |
| `notebooks/README.md` | what every `CHECK` establishes and why it is not circular. **Read this before changing either file** |
| `mc_config.m` | the shared geometry and protocol utilities: `pick_grid`, `rotate_about`, `load_protocol_file` |
| `mrtrix_io.m` | read and write MRtrix `.mif` / `.mih` images |
| `oct_path.m` | puts `SMI.m`, `helpers/` and the Octave shims on the path |
| `protocol/hcp_real_3shell.txt` | the acquisition in use, as tracked text |
| `stubs/` | Octave shims for `round(x,n)`, `discretize`, `datetime` |
| `check_manuscript_static.m` | static checks on the manuscript file: it parses, the scoring arrays are subscripted correctly, every `RUN{}` field read is written |
| `test_csd_arms.m` | the CSD arms alone in ~2 s, no `SMI.fit`. The regression test for the `-neg_lambda` bug |
| `measure_glyph_spread.m` | how much the drawn fODF glyph radius varies between noise realisations vs across SNR, ~4 min. Behind "README for Claude" section 6.6 |

The fODF machinery lives in `../helpers/`: `fODF_modulation_helpers.m` (forward
model, Watson, projection), `SMI_response_helpers.m` (kernel → zonal response,
glyphs) and `fODF_peak_score.m` (the one peak finder every arm goes through).

## Running it

```
cd notebooks
octave-cli --no-gui -q smi_manuscript_60deg.m     # SMOKE_TEST = true: minutes

cd ..
octave-cli --no-gui -q check_manuscript_static.m  # seconds, before any long run
octave-cli --no-gui -q test_csd_arms.m            # the MRtrix side alone, ~2 s
```

No Python and no data to download. MRtrix3 must be on the `PATH`. At its
manuscript settings the simulation is 42 `SMI.fit` calls and runs in hours;
`SMOKE_TEST = true` cuts it to one Lmax and three SNRs while still executing
every `CHECK`.

## The older campaign

`../Archive/deconv_pipeline/` holds the original Monte Carlo pipeline — the
Octave + MRtrix + Python arrangement that produced `Reports/deconv_tables.md`
and `Reports/REPORT_SMI_deconvolution_MonteCarlo.md`. It is archived rather than
deleted because those reports' numbers have no other provenance and
regenerating them on the real protocol is still an open task.

It reached 10,000 realisations and four crossing angles by splitting the work
across three languages and joining the arms by voxel index, which is exactly
what the manuscript file replaced. Prefer the manuscript file for anything new.

An even earlier version ran CSD and MSMT-CSD through dipy and reimplemented
`dwi2response dhollander`, `mrthreshold` and `amp2response` from the MRtrix
source. All of that is gone — MRtrix does it. See `git log` if it is ever
wanted back.

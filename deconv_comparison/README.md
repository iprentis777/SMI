# Constrained SMI vs SSST-CSD vs MSMT-CSD

The simulation behind the manuscript, in the design of Jeurissen et al. (2014):
synthesise noise-free signal vectors for crossing white matter fibres by forward
convolution, add complex Gaussian noise so the magnitude is Rician, and
deconvolve many independent realisations of each condition with each method.

**Everything here serves one file**, `notebooks/smi_wm_60deg.m`. It runs every
arm — SMI, SSST-CSD and both MSMT-CSD settings — on **one simulation, one noise
draw, one peak finder, in one scope**. There is no second forward model to keep
in step and no question about whether the arms saw the same data: they are the
same array.

Earlier simulations are in
[`../Archive/old_simulations/`](../Archive/old_simulations/), with the measured
defects that replaced them.

**CSD and MSMT-CSD are MRtrix3 3.0.4 itself.** `dwi2fod`, `dwiextract` and
`mrinfo` are the binaries, called as subprocesses. The only MRtrix behaviour
implemented locally is reading and writing its image format (`mrtrix_io.m`),
which `mrinfo` checks on every run.

**Everything here is simulation.** No result in this directory has touched a
patient scan.

## What is compared

| arm | what it is |
|---|---|
| SMI fixed | `SMI.get_plm_from_S_and_kernel` with a **fixed** kernel, the same one the CSD arms get as a response. No Rician bias correction, so the same noise model as MRtrix |
| SMI fitted | `SMI.fit`, estimating the kernel per voxel. Off by default — the expensive arm |
| SSST-CSD | `dwi2fod csd` on the top shell |
| MSMT def | `dwi2fod msmt_csd` at MRtrix's shipped `-neg_lambda 1e-10 -norm_lambda 1e-10` |
| MSMT tuned | the same, at `-neg_lambda 1 -norm_lambda 1e-3` |

**One response, built from one kernel, is handed unchanged to every arm.** That
is what makes the comparison answerable: the arms differ in their deconvolution
algorithm and in nothing else. The difference between the two SMI arms is
therefore the cost of estimating the kernel, measured rather than asserted.

Both MSMT settings run because the difference between them is **order
dependent** — a finding, not a bug to hide. `dwi2fod csd` and `dwi2fod msmt_csd`
do not ship comparable defaults (`-neg_lambda 1` against `1e-10`), so running
both "at their defaults" is not a like-for-like comparison and produced a wrong
result once already. See `Reports/REPORT_CSD_response_derivation.md`.

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
| `notebooks/smi_wm_60deg.m` | **the simulation.** Every arm, healthy WM, 60° crossing at 18 orientations, SNR swept, four figures |
| `notebooks/README.md` | what it asks, what the `CHECK` lines establish, and the conventions that silently ruin a result |
| `mc_config.m` | the shared geometry and protocol utilities: `pick_grid`, `rotate_about`, `load_protocol_file` |
| `mrtrix_io.m` | read and write MRtrix `.mif` / `.mih` images |
| `oct_path.m` | puts `SMI.m`, `helpers/` and the Octave shims on the path |
| `protocol/hcp_real_3shell.txt` | the acquisition in use, as tracked text |
| `stubs/` | Octave shims for `round(x,n)`, `discretize`, `datetime` |
| `test_csd_arms.m` | the CSD arms alone in ~2 s, no `SMI.fit`. The regression test for the `-neg_lambda` bug |

The fODF machinery lives in `../helpers/`: `fODF_sim_helpers.m` (forward
model, Watson, projection), `SMI_response_helpers.m` (kernel → zonal response,
glyphs) and `fODF_peak_score.m` (the one peak finder every arm goes through).

## Running it

```
cd notebooks
octave-cli --no-gui -q smi_wm_60deg.m      # HOURS at the shipped default

cd ..
octave-cli --no-gui -q test_csd_arms.m     # the MRtrix side alone, ~2 s
```

No Python and no data to download. MRtrix3 must be on the `PATH`.

⚠️ **The simulation ships at full size and runs for hours** — 18 orientations ×
50 reps × 16 SNR at Lmax 4, 6 and 8. Set `SMOKE_TEST = true` for a reduced run
in minutes, with every `CHECK` still executed and indicative numbers only. See
`notebooks/README.md` for why the full size is what the numbers were measured
at.

## The older campaign

`../Archive/deconv_pipeline/` holds the original Monte Carlo pipeline — the
Octave + MRtrix + Python arrangement that produced `Reports/deconv_tables.md`
and `Reports/REPORT_SMI_deconvolution_MonteCarlo.md`. It is archived rather than
deleted because those reports' numbers have no other provenance and
regenerating them on the real protocol is still an open task.

It reached 10,000 realisations and four crossing angles by splitting the work
across three languages and joining the arms by voxel index, which is exactly
what the single-file simulation replaced. Prefer `notebooks/smi_wm_60deg.m` for
anything new.

An even earlier version ran CSD and MSMT-CSD through dipy and reimplemented
`dwi2response dhollander`, `mrthreshold` and `amp2response` from the MRtrix
source. All of that is gone — MRtrix does it. See `git log` if it is ever
wanted back.

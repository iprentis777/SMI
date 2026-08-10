# The original Monte Carlo pipeline

Everything here produced `Reports/deconv_tables.md` and
`Reports/REPORT_SMI_deconvolution_MonteCarlo.md`. It is **archived, not
deleted**: those reports' numbers have no other provenance, and section 8 of
`README for Claude.md` still lists "regenerate the tables on the real HCP
protocol" as an open task, which would be a re-run of exactly these files.

**It is not the active comparison.** `deconv_comparison/notebooks/smi_manuscript_60deg.m`
replaced it and runs SMI, SSST-CSD and MSMT-CSD in one file on one noise draw.
Prefer that for anything new. This pipeline reached 10,000 realisations and four
crossing angles by splitting the work across Octave, MRtrix and Python, joining
the arms by voxel index across three languages rather than by being the same
array in one scope — which is the reason it was superseded.

## What is here

| file | what |
|---|---|
| `run_all.sh` | the whole campaign, from an empty `data/` |
| `setup_protocol.py` | the evaluation sphere; the only thing that needed `dipy` |
| `gen_phantom.m` | the synthetic voxel population the responses were estimated from |
| `gen_montecarlo.m` | the Monte Carlo signals and the SMI arm |
| `sweep_nonneg.m` | one kernel fit, six regularizer settings (report section 6) |
| `run_mrtrix.sh` | every MRtrix call: responses, fits, peaks, controls |
| `write_exact_response.py`, `kernel.py` | the exact kernel response in MRtrix format, and `K_l(b)` in Python |
| `binio.{m,py}` | the flat float64 Octave/Python exchange |
| `mrtrix_io.py` | the Python MRtrix image reader (the `.m` writer stayed in `deconv_comparison/`) |
| `score_mrtrix.py`, `score_sweep.py`, `tables.py` | scoring and the report tables |
| `figure_mc.py`, `figure_response.py` | `Figures/fodf_deconv_montecarlo.png` and `fodf_response_shview.png` |
| `dump_bases.m`, `check_mrtrix_basis.sh` | SMI's SH basis against MRtrix's, both `CS_phase` values |
| `mrtrix_responses/` | `dwi2response` outputs, estimated on the superseded synthetic protocol |
| `hcp_like_3shell.txt` | that superseded synthetic acquisition |
| `octave_test_stubs/` | graphics stubs from when Octave's `print` was thought unusable here. It is not — `figures_csd_arms`-style code writes real PNGs — so these are obsolete as well as archived |
| `smi_manuscript_60deg.ipynb` | the PR #22 three-arm Jupyter notebook, superseded by the `.m` |

## Running it again

The paths assume this directory sits beside `SMI.m` and `helpers/`, which it no
longer does. Anything here needs its `oct_path`/`addpath` lines pointed back at
the repository root before it will run. Nothing in the active manuscript path
imports any of it — that was checked before the move, and both test suites pass
without it.

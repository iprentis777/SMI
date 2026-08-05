# Monte Carlo comparison: constrained SMI vs MSMT-CSD vs SSST-CSD

The simulation behind `../Reports/REPORT_SMI_deconvolution_MonteCarlo.md`, in the design
of Jeurissen et al. (2014): synthesise noise-free signal vectors for crossing
white matter fibres by forward convolution, add complex Gaussian noise so the
magnitude is Rician, and deconvolve 10,000 independent realisations of each
condition with each method.

**CSD and MSMT-CSD are MRtrix3 3.0.4 itself** — `dwi2response`, `dwi2fod` and
`sh2peaks`, the binaries. The SMI fODF goes through the same `sh2peaks`, so peak
extraction is identical for every arm and is not something this package
implements. The only MRtrix behaviour implemented here is reading and writing
its image format (`mrtrix_io.m`, `mrtrix_io.py`), verified against `mrinfo` and
`mrconvert`.

**Everything here is simulation.** No result in this directory or in the report
has touched a patient scan — including the "real data" the response functions
are estimated from, which is a synthetic phantom. Section 2 of the report says
exactly what that does and does not buy.

## What is compared

| arm | what it is |
|---|---|
| SMI constrained | `SMI.fit`, `fODF_regularization.flag_nonneg = 1` at the shipped defaults, `CS_phase = 0` |
| SSST-CSD | `dwi2fod csd` on the b = 3 shell, `dwi2response tournier` response |
| MSMT-CSD | `dwi2fod msmt_csd` on all shells, `dwi2response dhollander` 3-tissue responses |

All three deconvolve at Lmax 6 and read the same image. The band limited ground
truth is written as an SH image too and run through the same `sh2peaks`, which
is where the "ceiling" column in the report comes from.

**`CS_phase = 0` matters.** At SMI's default of 1 the SH basis differs from
MRtrix's by `(-1)^m`, which is a 180 degree rotation about z of every fODF.
`check_mrtrix_basis.sh` measures that against MRtrix itself.

## Run order

`./run_all.sh` does all of it (`NREP=200 ./run_all.sh` for a quick version).
Step by step:

```
python3 setup_protocol.py                                   # gradient table
octave --eval "oct_path; gen_phantom(30,'p30')"             # response phantom
./run_mrtrix.sh responses                                   # dwi2response x3
./check_mrtrix_basis.sh                                     # SH basis vs MRtrix

octave --eval "oct_path; gen_montecarlo(30,10000,'snr30')"  # signal + SMI
./run_mrtrix.sh fit snr30                                   # dwi2fod + sh2peaks

python3 tables.py nf snr50:50 snr30:30 ...  | tee ../Reports/deconv_tables.md
python3 figure_mc.py ../Figures/fodf_deconv_montecarlo.png snr50:50 snr30:30 ...
python3 figure_response.py p30 ../Figures/fodf_response_shview.png

octave --eval "oct_path; sweep_nonneg(30,2000,'sw30')"      # report section 6
./run_mrtrix.sh sweep sw30 6 && python3 score_sweep.py sw30 30
./run_mrtrix.sh control snr50                               # report section 7
```

`nf` is the noise-free arm (`gen_montecarlo(1e4, 500, 'nf')`), used only to
normalise peak amplitudes so methods whose fODFs live on different scales can
still be compared.

## Files

| file | what |
|---|---|
| `setup_protocol.py` | the HCP-like 3-shell gradient table (the only thing that uses dipy) |
| `mc_config.m` | **every constant of the experiment, in one place.** `gen_montecarlo.m` and `sweep_nonneg.m` both read from it, so the main run and the sweep cannot disagree about what is being simulated. Start here to see what the experiment *is* |
| `gen_phantom.m` | the synthetic voxel population the responses are estimated from |
| `gen_montecarlo.m` | the Monte Carlo signals and the SMI arm |
| `sweep_nonneg.m` | one kernel fit, six regularizer settings (report section 6) |
| `run_mrtrix.sh` | every MRtrix call: responses, fits, peaks, controls |
| `write_exact_response.py` | the exact kernel response in MRtrix response format |
| `mrtrix_io.{m,py}` | read and write MRtrix `.mif` / `.mih` images |
| `score_mrtrix.py` | matches `sh2peaks` output to the true axes; the ACC |
| `score_sweep.py`, `tables.py`, `figure_*.py` | tables and figures |
| `check_mrtrix_basis.sh` | SMI's SH basis against MRtrix's, both `CS_phase` values |
| `dump_bases.m` | writes SMI's SH basis for that check |
| `binio.{m,py}`, `kernel.py` | Octave/Python exchange, and the SM kernel in Python |
| `stubs/` | Octave shims for `round(x,n)`, `discretize`, `datetime` |
| `octave_test_stubs/` | graphics stubs, so the plotting examples run headless |

`binio.{m,py}` and `kernel.py` are copies of the files of the same name in
`freewater_comparison/` on the unmerged `claude/freewater-simulations` branch.
If that branch lands they should be de-duplicated; they are copied rather than
referenced so this package does not depend on an unmerged branch.

Requires **MRtrix3 3.0.4**, Octave (`statistics`) and Python
(`numpy matplotlib`, plus `dipy` for `setup_protocol.py` alone). Intermediate
arrays land in `data/` and `mrtrix/`, neither of which is tracked.

## Planned extension: edematous kernel

A future simulation condition should consider a kernel in an edematous
environment. This is a scope marker only: no edema model, parameterization,
response-estimation strategy, or result is defined yet.

## Superseded

An earlier version of this package ran CSD and MSMT-CSD through dipy and
reimplemented `dwi2response dhollander`, `mrthreshold` and `amp2response` from
the MRtrix source. That is all removed — MRtrix does it. Its conclusions about
MSMT-CSD reproduced under the real thing; see `git log` if the reimplementation
is ever wanted back.

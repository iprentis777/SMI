# Monte Carlo comparison: constrained SMI vs MSMT-CSD vs SSST-CSD

The simulation behind `../REPORT_SMI_deconvolution_MonteCarlo.md`, in the design
of Jeurissen et al. (2014): synthesise noise-free signal vectors for crossing
white matter fibres by forward convolution, add complex Gaussian noise so the
magnitude is Rician, and deconvolve 10,000 independent realisations of each
condition with each method.

**Everything here is simulation.** No result in this directory or in the report
has touched a patient scan. The response functions are estimated by a
reimplementation of `dwi2response dhollander` from a synthetic phantom, not from
real data — see section 2 of the report for exactly what that does and does not
buy.

## What is compared

| arm | what it is |
|---|---|
| SMI constrained | `SMI.fit` with `fODF_regularization.flag_nonneg = 1` |
| SMI unconstrained | the same kernel, deconvolved by plain LLS |
| SSST-CSD | dipy `ConstrainedSphericalDeconvModel`, b = 3 shell, `dwi2response tournier` response |
| MSMT-CSD | dipy `MultiShellDeconvModel`, all shells, `dwi2response dhollander` responses |

The two SMI arms share one kernel fit, so the only difference between them is
the deconvolution. All four are deconvolved at Lmax 6, read the same DWI bytes,
and are scored by one peak finder on one sphere.

## Run order

```
python3 setup_protocol.py                     # protocol + evaluation sphere
octave --eval "oct_path; gen_phantom(30,'p30')"          # response phantom
python3 dhollander.py p30                     # dhollander / tournier / fa responses

# per SNR arm (these are independent and can run concurrently)
octave --eval "oct_path; gen_montecarlo(30,10000,'snr30')"
python3 run_csd.py snr30 p30

python3 check_conventions.py snr30            # basis + convolution checks
python3 score.py snr30 30 nf                  # tables
python3 figure_mc.py ../fodf_deconv_montecarlo.png snr50:50 snr30:30 ...
python3 figure_response.py p30 ../fodf_response_shview.png
```

`nf` is the noise-free arm (`gen_montecarlo(1e4, 500, 'nf')`), used only to
normalise peak amplitudes so that methods whose fODFs live on different scales
can still be compared.

## Files

| file | what |
|---|---|
| `setup_protocol.py` | HCP-like 3-shell protocol and the shared evaluation sphere |
| `gen_phantom.m` | the synthetic voxel population the responses are estimated from |
| `mrtrix_ops.py` | `mrthreshold`, `mrthreshold -top`, `amp2response` |
| `mtcsd.py` | 2-tissue MSMT used only inside the Dhollander single-fibre step |
| `dhollander.py` | `dwi2response dhollander`, plus `tournier` and `fa` |
| `gen_montecarlo.m` | the Monte Carlo signals and both SMI arms |
| `run_csd.py` | SSST-CSD and MSMT-CSD through dipy |
| `peaks.py` | the one peak finder, with sub-vertex refinement |
| `score.py` | the tables |
| `check_conventions.py` | SH basis and convolution checks |
| `figure_*.py` | the figures |
| `binio.{m,py}`, `kernel.py` | exchange layer and the SM kernel in Python |
| `stubs/` | Octave shims for `round(x,n)`, `discretize`, `datetime` |
| `octave_test_stubs/` | graphics stubs, so the plotting examples can be run headless |

`binio.{m,py}` and `kernel.py` are copies of the files of the same name in
`freewater_comparison/` on the unmerged `claude/freewater-simulations` branch.
If that branch lands they should be de-duplicated; they are copied rather than
referenced so this package does not depend on an unmerged branch.

Requires Octave (`statistics`) and Python (`numpy scipy dipy cvxpy matplotlib`).
`cvxpy` is not optional: dipy's MSMT-CSD QP and `amp2response`'s constraints both
need it. Intermediate arrays land in `data/`, which is not tracked.

## Extra arms

Beyond the four-method comparison:

| script | what |
|---|---|
| `sweep_nonneg.m` + `score_sweep.py` | the `lambda_nonneg` sweep of report section 6 — one `SMI.fit`, six deconvolutions from the same kernel |
| `control_exact.py` + `score_control.py` | CSD and MSMT-CSD given the exact kernel response instead of an estimated one (report section 7) |
| `dump_bases.m` | SMI's SH basis with and without the Condon-Shortley phase, for the MRtrix basis check in `check_conventions.py` |
| `tables.py` | emits the report's result tables as markdown, so no number in the report is transcribed by hand |

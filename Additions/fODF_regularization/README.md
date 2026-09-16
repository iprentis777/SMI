# Addition: regularized fODF deconvolution

**Status:** active, opt-in, off by default.
**Option:** `options.fODF_regularization`
**Stock NYU SMI:** has none of this. The deconvolution there is an unregularized
least-squares fit, which is what you still get unless you set `flag_nonneg = 1`.

## The problem

Once the kernel is known, the fODF coefficients `plm` are obtained by
deconvolving the kernel from the signal. That is ill-conditioned: the kernel
rotational invariants `K_l` decay quickly with `l`, so the high-order `plm` are
dominated by noise and the estimated fODF develops large negative lobes.

Negative fODF amplitude is physically meaningless. On the synthetic two-fibre
benchmark at SNR 30, three shells, Lmax 6:

| | relative fODF error | negative mass |
|---|---|---|
| unregularized | 0.615 | 0.159 |
| non-negativity constrained | **0.199** | **0.031** |
| *ground truth, for reference* | — | 0.0034 |

## What it does

A non-negativity constraint, following constrained spherical deconvolution
([Tournier et al., NeuroImage 2007](https://doi.org/10.1016/j.neuroimage.2007.02.016)).
The fODF is first estimated with a low-order unconstrained fit, then refined:
at each iteration the directions where the fODF falls below `tau*mean(fODF)` are
collected, a penalty on their amplitude is added to the least-squares problem,
and this repeats until that set of directions stops changing — typically 4–6
iterations.

```matlab
options.flag_fit_fODF = 1;
options.fODF_regularization.flag_nonneg   = 1;    % default 0 (off)
options.fODF_regularization.lambda_nonneg = 1;    % default 1
options.fODF_regularization.tau           = 0.1;  % default 0.1
options.fODF_regularization.Ndirs         = 300;  % default 300
options.fODF_regularization.Niter         = 50;   % default 50 (max)
options.fODF_regularization.Lmax_init     = 4;    % default 4

out = SMI.fit(dwi, options);

out.fODF_regularization.Niterations     % iterations used, per voxel
out.fODF_regularization.Nnegative_dirs  % constrained directions, per voxel
out.fODF_regularization.flag_converged  % per voxel
```

`lambda_nonneg` is dimensionless: the regularization block is rescaled by the
norm of the rows of each voxel's design matrix, so one value means the same
thing across voxels and protocols.

## The one setting worth arguing about

**`lambda_nonneg` defaults to 1, and two measurements disagree about that.**

It was briefly changed to 10 on the strength of a sweep minimizing the relative
L2 error of the fODF over the sphere, then **changed back**, because a much
larger Monte Carlo scoring peak orientation found:

| `lambda_nonneg` | 45° crossings resolved | angular correlation (1 fibre / 15° / 45° / 60°) |
|---|---|---|
| **1** | **55%** | 0.980 / 0.986 / 0.966 / 0.972 |
| 3 | 0.2% | — |
| 10 | **0.0%** at every SNR, including noise-free | 0.930 / 0.947 / 0.887 / 0.903 |

The two scores measure different things. L2 error over the sphere is dominated
by the isotropic part and by negative mass, so it happily trades angular
resolution for smoothness. Peak orientation and fibre count are what a
tractography algorithm actually consumes.

Spurious peaks are **already fully suppressed at 1** — every constrained setting
sits at exactly 0.000 spurious peaks per voxel, while turning the constraint off
puts 0.026 per voxel into a 45° crossing. Since nothing above 1 buys further
suppression and everything above 1 costs angular resolution, **1 is the
default.** Raise it toward 10 only if smoothness of the whole fODF matters more
to you than resolving crossings.

## Does it flatten the fODF?

No — and this is the counterintuitive part. The **unregularized** fODF is not
the reference height: its peaks are **inflated by 12%**, because the peak is a
maximum over directions and noise biases a maximum upward.

| | peak(estimate) / peak(truth) |
|---|---|
| unregularized | 1.116 |
| `lambda_nonneg = 1` | 1.019 |
| `lambda_nonneg = 10` | 1.005 |

The constraint **removes the noise-driven inflation** rather than shrinking the
fODF below truth. (Only when pushed far past its optimum, at 100, does the ratio
rise again to 1.30 — a constraint that strong forces the fODF to zero over most
of the sphere and concentrates the remaining mass into narrow spikes.)

## Known limitation

The constraint is what suppresses high-`l` power in the recovered fODF. Together
with error in the estimated kernel — whose `K_l` at high `l` is small and very
sensitive — it is why `pl4` comes back at roughly a third of its true value
(median 0.26 in single-fibre WM against a true 0.71). This was for a long time
blamed on Tikhonov damping; that was measured wrong, and is part of why the
damping was removed. If you need faithful high-`l` power, this is the term
responsible.

## Where to look next

- Full documentation and all result tables: root [`README.md`](../../README.md),
  "Regularized fODF deconvolution".
- Runnable example, no data needed:
  [`examples/example_fODF_regularization.m`](../../examples/example_fODF_regularization.m)
- Parameter sweep:
  [`examples/example_fODF_regularization_sweep.m`](../../examples/example_fODF_regularization_sweep.m)
- Measurement report:
  [`Reports/REPORT_fODF_regularization_sweep.md`](../../Reports/REPORT_fODF_regularization_sweep.md)
- The retired sibling:
  [`Archive/patch_history/fODF_tikhonov/`](../../Archive/patch_history/fODF_tikhonov/)

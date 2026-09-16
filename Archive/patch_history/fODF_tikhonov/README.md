# Retired: Tikhonov damping of the fODF deconvolution

**Status:** retired. Removed from the toolbox; no code path reads it.

## The question it was meant to answer

The fODF deconvolution is ill-conditioned: the kernel rotational invariants
`K_l(b)` decay quickly with `l`, so the high-order `plm` are dominated by noise
and the resulting fODF has large negative lobes at clinical SNR.

Tikhonov damping added a penalty

```text
lambda_tikhonov^2 * || Gamma * plm ||^2
```

to the least-squares problem, with `Gamma` either the identity or the
Laplace-Beltrami matrix `diag(l(l+1))` normalised by its maximum. The reasoning
was that it would suppress exactly the coefficients the kernel attenuates most.
It shipped alongside the non-negativity constraint, was disabled by default
(`lambda_tikhonov = 0`), and was introduced in
`Patches/0001fODFregularization.patch`.

## Why it is retired: inert where it was used, harmful where it was not

The reasoning above was never confirmed by a measurement. Three measurements
already in this repository contradict it, and together they leave no weight at
which the term earns its place.

**At the weights anyone actually used (`lambda <= 0.8`, default 0) it is
inert:**

1. **It does not touch the bands it was supposed to touch.** Swept from 0 to
   0.8 at three noise levels, the high-`l` bands were identical to three
   decimal places. Recorded in `SMI.m`, in the comment above
   `fODF_ModulationWeight` that now reads "an earlier version of this comment
   blamed Tikhonov damping for the lost `l=4` power".
2. **It does not move an angular result.** In the 60 degree crossing
   simulation, `lambda_tikhonov = 0.3` against `0` moved a 45 degree error from
   **21.28 to 21.27 degrees**.

**At larger weights it is not inert — it is harmful:**

3. **It is the thing that flattens the fODF**, strongly and monotonically. From
   the sweep's peak amplitude ratio `peak(estimate)/peak(truth)`: at
   `lambda_tikhonov = 3` the peaks retain **55%** of their true height, and at
   `10` only **27%**. That is the mechanism behind the large errors in the high
   weight columns of the sweep's accuracy table. The non-negativity constraint,
   by contrast, does not flatten the fODF — it removes a noise-driven *inflation*
   of about 12%.

So the honest summary is not simply "it did nothing". It did nothing in the
range it shipped in, and everything it did outside that range made the result
worse. There is no weight at which it was measured to help.

The first measurement was made while investigating a *different* question — why
`l = 4` power was being lost — and it overturned the then-current explanation.
The cause of that loss is the **non-negativity constraint together with error in
the estimated kernel**, whose `K_l` at high `l` is small and very sensitive.
That conclusion is unchanged by this removal; the removal confirms it.

## What was learned

- The high-order power loss in SMI fODFs is a **non-negativity** effect, not a
  damping effect. Anyone re-investigating it should start at `lambda_nonneg`.
- Non-negativity is the only fODF regularizer in SMI that demonstrably improves
  a result, which makes it the only one worth tuning. It is also the term with
  real and measurable side effects, so "less regularization" is not free.
- **Peak height and accuracy come apart.** The unregularized fODF's peaks are
  inflated 12% by noise, so "taller" is not "more faithful". Damping shrank
  peaks below the truth; the constraint removed the inflation without
  overshooting. Any future regularizer should be scored on both.
- Removing this leaves SMI **less** regularized than the CSD arms it is
  compared against, not equally: `dwi2fod csd` ships `-norm_lambda 1`, a
  penalty on the norm of the solution, which is Tikhonov with `Gamma = I`.
  Any comparison against CSD should state this rather than imply parity.

## Why it is not recommended pipeline guidance

It is not available at all. An unmotivated, measured-inert knob is a question a
reviewer will ask for no measured benefit.

**A caller that still sets `lambda_tikhonov` is silently ignored.**
`SMI.fODF_RegularizationDefaults` fills defaults for the fields it knows and
does not reject unknown ones, so an old script will run without error and
without damping. If you are re-running an old analysis and comparing against
old numbers, that is the thing to check first.

## Provenance

| file | what |
|---|---|
| `0001-remove-tikhonov-damping-solver.patch` | the removal from `SMI.m`, `SMI_freeL0.m` and `helpers/fODF_free_l0_deconv.m` — the implementation itself |
| `0002-remove-tikhonov-damping-callers.patch` | the removal from the examples, tests and simulation notebooks that set the option |
| `../../../Patches/0001fODFregularization.patch` | where it was introduced, with its original rationale |
| `../../../Reports/REPORT_fODF_regularization_sweep.md` | the sweep that measured the regularization weights |

Both patches are plain `git diff` output. To restore the feature, apply them in
reverse from the repository root:

```sh
git apply -R Archive/patch_history/fODF_tikhonov/0002-remove-tikhonov-damping-callers.patch
git apply -R Archive/patch_history/fODF_tikhonov/0001-remove-tikhonov-damping-solver.patch
```

They were generated against the commit that removed the feature, so a reverse
apply onto a later tree may need `--3way`.

## Before restoring it

Restoring this is only worth doing with a measurement that the two above did
not make. Both were made at a fixed, known kernel on healthy white matter. A
claim that damping helps would need to show an effect that survives at a
*fitted* kernel, or on a tissue where `K_l` behaves differently — and it would
need to beat simply moving `lambda_nonneg`, which is the knob that was measured
to do something.

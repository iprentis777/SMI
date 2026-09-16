# Experimental: free `l = 0` deconvolution

**Status:** experimental, archived. Never part of the active toolbox.

Everything here is additive — it was deliberately built as **separate files**
rather than as changes to `SMI.m`, so `SMI.m` never depended on it and removing
it from the active tree changes nothing about a normal fit.

## The question

SMI's deconvolution imposes `p_00 = 1`, so the fODF integrates to 1 in every
voxel. Two consequences, both of which matter for tractography:

- **the fODF carries no density information at all** — a CSF voxel and a
  coherent white matter voxel have equal mass;
- its isotropic floor is the constant `1/(4*pi) = 0.0796`, which sits **above**
  MRtrix's default `iFOD2 -cutoff 0.05`. An SMI fODF therefore passes the
  tractography termination test in *every voxel of the brain, CSF included*.
  Measured: SMI leaves CSF at 0.32 where MSMT-CSD leaves it at 0.028.

This work asked what happens if that constraint is lifted and every coefficient
**including `l = 0`** is estimated from the data.

It is the same problem that
[anisotropy modulation](../patch_history/fODF_modulation/) attacked from the
other direction — modulation kept `p_00 = 1` and reweighted the result, while
this estimates `p_00` directly. Neither shipped. The problem is real and still
open.

## What is here

| file | what |
|---|---|
| `SMI_freeL0.m` | a full fork of `SMI.m` carrying the `free_l0` option and a one-step active set on `p_00` with `p00_bounds` |
| `fODF_free_l0_deconv.m` | **both conventions in one solver**, so the only difference between the two arms is the convention itself |
| `smi_free_l0_experiment.m` | the comparison, with a check that `'fixed'` mode reproduces `SMI.get_plm_from_S_and_kernel` to machine precision |

The single-solver design is the good idea worth keeping: if `'fixed'` mode does
not reproduce the shipped deconvolution exactly, nothing the experiment reports
is worth reading, and Step 2 checks exactly that before anything else runs.

## Why it is archived

Not because it was measured wrong — it was never measured to a conclusion. It
is archived because:

- `SMI_freeL0.m` is a **whole-file fork of `SMI.m`** (both over a megabyte), and
  a fork drifts. Two copies of a solver is how a fix lands in one and not the
  other;
- its only driver was one experiment script, so nothing in the repository
  exercised it;
- it was explicitly landed as "EXPERIMENTAL: optional free `l = 0`
  deconvolution, **in additive files only**" and never graduated.

Keeping a second megabyte-scale solver in the active tree implies a support
commitment the work had not earned.

## If you pick this up again

`SMI_freeL0.m` is a fork of `SMI.m` **as it was when the fork was taken**, and
`SMI.m` has moved since. Two features were removed from both at once, so the
fork is current on those:

- Tikhonov damping — removed from `SMI_freeL0.m` and
  `fODF_free_l0_deconv.m` as well. See [`../patch_history/fODF_tikhonov/`](../patch_history/fODF_tikhonov/).
- anisotropy modulation — likewise. See
  [`../patch_history/fODF_modulation/`](../patch_history/fODF_modulation/).

Anything else that changed in `SMI.m` after the fork point is **not** reflected
here. Diff the two before trusting this file:

```sh
diff <(sed 's/[[:space:]]*$//' SMI.m) \
     <(sed 's/[[:space:]]*$//' Archive/free_l0/SMI_freeL0.m)
```

The right way to revive this is almost certainly **not** to resurrect the fork.
It is to take the single-solver approach of `fODF_free_l0_deconv.m` — both
conventions in one function, with an exactness check against the shipped path —
and put the option into `SMI.m` itself.

Note also that the experiment script calls `fODF_free_l0_deconv()` by bare name
and expects it on the path. From this archived location that needs an
`addpath` to this directory, not to `helpers/`.

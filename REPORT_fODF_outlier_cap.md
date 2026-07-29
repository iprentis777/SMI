# Post hoc outlier cap on the fODF

Measurement report for `SMI.fODF_OutlierDefaults`, `SMI.cap_fODF_outliers` and
`SMI.neighbour_median`. Everything below is measured; the test harnesses are
`test_fODF_outlier_cap.m` (the method in isolation) and `test_SMI_outlier_cap.m`
(the flag on `SMI.fit`), both of which are in this PR and need no data.

The numbers used to choose the ceiling in section 4, and the peak distributions
in section 5, come from a free-water simulation that is deliberately NOT part of
this PR — it lives on `claude/freewater-simulations` and will land separately.
The cap itself has no dependency on it: nothing in `SMI.m` or in either test
harness references that code, so this PR stands alone.

---

## 1. What it does

A voxel is flagged when **either**

```
peak > 10^orders * (median peak of its 26 in-mask neighbours)     % relative
peak > ceiling                                                    % absolute
```

and is then rescaled **down** to `min(neighbourhood median, ceiling)`. The
operation is strictly one sided: a voxel is never raised. It can remove
spurious amplitude but can never invent fibre density where there was none,
which is the property that makes it safe to leave on.

Defaults: `orders = 1` (one order of magnitude), `ceiling = 1`,
`min_neighbours = 6`, `connectivity = 26`, `Ndirs = 500`.

---

## 2. Why this is not a tissue type criterion

Every per voxel weight in `REPORT_fODF_modulation.md` was rejected, or nearly
rejected, because it keyed on something edema shares with CSF. The
neighbourhood test cannot: it is **relative**, and edema is spatially
**contiguous**. An edematous voxel's neighbours are edematous too, so the local
median moves with it and the ratio does not. A whole region being uniformly
bright or dim never trips this; only isolated spikes do.

That is a structural property of the statistic, not a threshold that happens to
be tuned well, and it is the first candidate in this project that gets the edema
constraint for free. It is verified directly: a contiguous 3x3x3 block raised by
a uniform factor is **never** flagged, at any of the settings tested, while
isolated spikes in the same volume always are.

---

## 3. Three things that are easy to get wrong

**The scale factor is not `target/peak`.** `out.plm` is stored with `p_00 = 1`,
so the isotropic floor is a fixed `1/(4*pi) = 0.0796` that is *not* part of the
blow-up: the excess lives entirely in `l >= 2`. Bringing the peak to `T` means
scaling only the anisotropic part,

```
s = (T - 1/(4*pi)) / (peak - 1/(4*pi))
```

Using `T/peak` would shrink the floor as well and undershoot. This is the same
trap flagged for fODF rescaling in `README for Claude.md` section 3.

**Median, not mean and standard deviation.** The median has a 50% breakdown
point; a mean and SD have none. A single 1e13 voxel among 26 neighbours inflates
the SD to ~2e12, so the offender sits 0.2 SD above the mean and is never
flagged. Worse, ventricles are contiguous, so several neighbours are usually
blown up together and a mean/SD rule is blind by construction.

**Scale the whole `l >= 2` block by one scalar.** That preserves peak
orientation *exactly*: the capped peak lands on the identical sphere vertex for
spikes of 40x, 5e3x, 1e11x and 1e12x. Clipping coefficients independently would
not, and can swing peak orientations.

---

## 4. On the ceiling, and why it is not 0.4-0.5

Measured at Lmax 6 on the free-water simulation branch (`claude/freewater-simulations`, not part of this PR):

| | peak |
|---|---|
| hard physical max of a band-limited fODF, `sum_{l even <= 6} (2l+1)/(4*pi)` | **2.228** |
| ground truth single fibre, Watson kappa 16 | 1.459 |
| what SMI actually recovers for single-fibre WM (median) | 0.84-0.86 |
| **fraction of simulated WM voxels above 0.5** | **60-64%** |
| fraction of simulated WM voxels above 1.0 | **0%** |

A ceiling of 0.4-0.5 sits in the **middle** of the white matter distribution and
would flatten most of it. The default of 1 clears every legitimate voxel with
headroom.

**But 1 is an empirical ceiling, not a physical one.** It is *below* the true
single-fibre peak of 1.459, and is only safe because SMI under-recovers the high
`l` bands — measured noise-free, `p6` comes back at 28% of truth and `p4` at
62%, so reconstructed peaks top out around 0.88. If the deconvolution is ever
sharpened (for instance by lowering `lambda_nonneg`, which the regularization
sweep on that branch shows recovers 45 degree crossings), peaks will rise and
this ceiling will start clipping real voxels. `2.228` is the
value that is a genuine bound.

This was not a theoretical worry. The first run of the test harness used ground
truth Watson kappa 16 fODFs, peak 1.445, and **all 720 voxels tripped the
ceiling of 1**. The harness now applies SMI's measured per band under-recovery
so the phantom matches what the pipeline actually produces.

---

## 5. Does it have anything to correct?

Peak fODF amplitude, SNR 15, Lmax 6, from the free-water simulation branch (`claude/freewater-simulations`, not part of this PR):

| regularization | CSF median | CSF p99 | CSF max | WM median | max / WM med |
|---|---|---|---|---|---|
| **nonneg 10 + tikhonov 0.3 (shipped)** | 0.395 | 0.938 | **0.977** | 0.561 | **1.7x** |
| nonneg 10 + tikhonov 0 | 0.392 | 0.940 | 0.978 | 0.563 | 1.7x |
| nonneg OFF + tikhonov 0.3 | 1.873 | 9.885 | 10.091 | 0.919 | 11.0x |

**With the shipped regularization there is no outlier population in simulation
at all** — the worst CSF voxel in 40 realisations is 1.7x the median WM peak.
Blow-ups only appear once the non-negativity constraint is off.

So on simulated data this cap is insurance, not a fix, and `Ncap = 0` is the
expected result. It is worth running on real data precisely because a non-zero
count would mean the regularization is not doing there what it does here — and
that is a more important finding than the cap itself. It also answers the open
question in `README for Claude.md` section 11 about whether the stray giant
glyphs are inside `mask3D`: `out.fODF_outlier.flagged` says exactly where they
are.

---

## 6. Verification

`test_fODF_outlier_cap.m` — the method on a 12x12x5 phantom of realistic
SMI-amplitude fODFs, with planted spikes at 50x, 1e4x and 1e12x, a spike inside
a near-isotropic (CSF analogue) region, and a contiguous bright block:

| check | result |
|---|---|
| all 4 planted spikes flagged | PASS |
| contiguous bright block NOT flagged | PASS (0 of 27) |
| only flagged voxels changed | PASS (4 changed, 4 flagged) |
| capped peak lands on `min(neighbourhood median, ceiling)` | PASS |
| capped peak lands on the identical sphere vertex | PASS (0.000 deg) |
| nothing left above the ceiling | PASS |
| every unflagged voxel bit identical | PASS (max abs diff 0 over 716) |
| both rules verified in isolation (ceiling disabled) | PASS |

Orientation is asserted by **sphere vertex identity**, not by an angle
tolerance: `acosd` cannot resolve below about 1e-5 degrees near a dot product of
1 and returns floating point noise there, so an angle threshold would be flaky.
The residual 1.2e-06 degrees reported below is `acosd(1-eps)` on the *same*
vertex, not a real rotation.

`test_SMI_outlier_cap.m` — the flag on a real `SMI.fit`, 8x8x4 at SNR 30:

| check | result |
|---|---|
| flag off is bit identical to never setting the option (plm, pl, kernel) | PASS |
| flag off leaves `out.plm_capped` absent | PASS |
| flag on leaves `out.plm` identical to flag off | PASS |
| 1e11x spike capped, peak on the identical sphere vertex | PASS (1.2e-06 deg, see below) |
| contiguous bright block not capped | PASS |
| modulation applied to the CAPPED fODF | PASS (max abs diff 0) |

---

## 7. Usage

```matlab
options.flag_fit_fODF = 1;
options.fODF_outlier.flag_cap = 1;     % default 0, off
out = SMI.fit(dwi, options);

out.plm_capped                          % corrected coefficients
out.fODF_outlier.flagged                % which voxels were capped
out.fODF_outlier.scale                  % the scalar applied, 1 elsewhere
out.fODF_outlier.peak_before            % peak amplitude map before the cap
out.fODF_outlier.Ncap                   % counts, also written to the log
```

or post hoc on an already fitted `out`:

```matlab
[plm_capped, info] = SMI.cap_fODF_outliers(out);
[plm_capped, info] = SMI.cap_fODF_outliers(out, struct('ceiling',2.228,'orders',1.5));
```

`out.plm`, `out.pl` and `out.kernel` are **identical whether the flag is on or
off**, matching the convention the regularization and modulation already follow,
so a capped and an uncapped run can be compared directly.

**Order of operations.** If both this and the modulation are enabled, the cap
runs first and the modulation is applied to the corrected fODF. That ordering is
forced: the cap works on absolute amplitudes and modulation rescales them, the
same reason peak truncation has to precede modulation.

---

## 8. Exactly what changed in `SMI.m`

| lines | what |
|---|---|
| 155-187 | header documentation of the `options.fODF_outlier` block |
| 340-351 | option parsing, including the guard that rejects `flag_cap = 1` with `flag_fit_fODF = 0` |
| 613-637 | the call site in `SMI.fit`, before the modulation, writing `out.plm_capped` and `out.fODF_outlier` |
| 745-763 | the log file entries, including the explicit `none` line when the flag is off |
| 1433-1477 | `SMI.fODF_OutlierDefaults` — option defaults and validation |
| 1478-1668 | `SMI.cap_fODF_outliers` — the cap itself |
| 1669-1701 | `SMI.neighbour_median` — masked neighbourhood median, NaN padded |

Line numbers are as of the commit that introduced this.

**One existing line was modified**, which is a departure from the modulation
work where every hunk was purely additive:

```
-   [sh_mod,w_mod,info_mod] = SMI.modulate_fODF(out,fODF_modulation);
+   [sh_mod,w_mod,info_mod] = SMI.modulate_fODF(out_mod,fODF_modulation);
```

`out_mod` is `out` with `plm` replaced by `plm_capped` when the cap ran, and is
`out` itself otherwise, so **a fit that does not enable the cap is unchanged**.
The change is what makes the two features compose in the documented order.

`SMI.neighbour_median` takes the median by sorting rather than with an
`'omitnan'` flag. That needs no toolbox and no version-dependent nanflag
support, and it is why the whole thing runs under Octave as well as MATLAB.

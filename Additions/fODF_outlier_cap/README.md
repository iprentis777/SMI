# Addition: post hoc fODF outlier cap

**Status:** active, opt-in, off by default.
**Option:** `options.fODF_outlier`
**Stock NYU SMI:** has none of this.

## The problem

Occasionally a single voxel returns an fODF amplitude orders of magnitude above
everything around it. The deconvolution gain `g_2 = 1/||K_2||` blows up where the
kernel is nearly isotropic, so a near-CSF voxel can reach enormous amplitudes —
unregularized at SNR 15 the simulation reaches **2.2e14**. Those glyphs dominate
any display and any amplitude threshold downstream.

## What it does

A voxel is flagged when **either**

```text
peak > 10^orders  x  (median peak of its in-mask neighbours)
peak > ceiling
```

and is then scaled **down** to `min(neighbourhood median, ceiling)`.

Strictly one-sided: a voxel is never raised, so this can remove spurious
amplitude but can never invent fibre density where there was none.

```matlab
options.flag_fit_fODF = 1;
options.fODF_outlier.flag_cap       = 1;    % default 0 (off)
options.fODF_outlier.orders         = 1;    % default: one order of magnitude
options.fODF_outlier.ceiling        = 1;    % default: absolute peak ceiling
options.fODF_outlier.min_neighbours = 6;    % default
options.fODF_outlier.connectivity   = 26;   % default, or 6

out = SMI.fit(dwi, options);

out.plm_capped                % corrected coefficients, same convention as out.plm
out.fODF_outlier.flagged      % which voxels were capped
out.fODF_outlier.scale        % scalar applied per voxel, 1 elsewhere
out.fODF_outlier.peak_before  % peak amplitude map before the cap
```

`out.plm`, `out.pl` and `out.kernel` are **identical whether the flag is on or
off** (verified, difference exactly 0). Also available post hoc as
`SMI.cap_fODF_outliers(out, options)`.

## Why this is not a tissue-type criterion

This matters for edema, where tissue-type thresholds fail.

The neighbourhood test is **relative**, and edema is spatially **contiguous** —
an edematous voxel's neighbours are edematous too, so the local median moves with
it and the ratio does not. A whole region being uniformly bright or dim never
trips this; only isolated spikes do.

That is a structural property of the statistic, not a well-tuned threshold, and
it is verified directly: a contiguous block raised by a uniform factor is never
flagged, while isolated spikes in the same volume always are.

## Three details that are easy to get wrong

1. **The scale factor is not `target/peak`.** `plm` is stored with `p_00 = 1`, so
   the isotropic floor `1/(4*pi) = 0.0796` is not part of the blow-up — the
   excess is entirely in `l >= 2`. The correct factor is
   `s = (T - 1/(4*pi)) / (peak - 1/(4*pi))`, applied to the `l >= 2` block.
   `T/peak` shrinks the floor too and undershoots.
2. **Median, not mean and standard deviation.** The median has a 50% breakdown
   point; an SD has none. One `1e13` voxel among 26 neighbours inflates the SD
   until the offender sits 0.2 SD above the mean and is never flagged — and
   blown-up voxels usually arrive in contiguous clusters.
3. **Scale the whole `l >= 2` block by one scalar.** That preserves peak
   orientation *exactly* — measured 0.000 deg for spikes of 40x, 5e3x and 1e11x.
   Clipping coefficients independently does not, and can swing orientations.

## Choosing the ceiling

| | peak |
|---|---|
| hard physical max of a band-limited fODF at Lmax 6, `sum (2l+1)/(4*pi)` | **2.228** |
| ground truth single fibre, Watson `kappa = 16` | 1.459 |
| what SMI actually recovers for single-fibre WM | 0.84–0.86 |
| fraction of simulated WM voxels above 0.5 | **60–64%** |
| fraction of simulated WM voxels above 1.0 | 0% |

A ceiling of 0.4–0.5 sits in the **middle** of the white matter distribution and
would flatten most of it. The default `1` clears every legitimate voxel with
headroom — but it is an **empirical** ceiling, not a physical one. It sits below
the true single-fibre peak of 1.459 and is only safe because SMI under-recovers
the high `l` bands (`p6` returns at 28% of truth). **If the deconvolution is ever
sharpened, raise it.** `2.228` is the only value that is a genuine bound.

## Expect zero caps on a regularized fit

With the regularization enabled the simulation has no outlier population at all:
at SNR 15 the worst CSF voxel is 1.7x the median WM peak. Blow-ups only appear
once the non-negativity constraint is off (max/median 11.0x).

So **`Ncap = 0` is the expected result**, and a non-zero count on real data is
itself the finding — it means the regularization is not behaving there as it does
in simulation. Treat the cap as a diagnostic first and a correction second, and
look at *where* the flagged voxels are: in the ventricles means the deconvolution
is blowing up there; on the brain edge means the mask erosion is not aggressive
enough.

## Order of operations

The cap works on **absolute amplitudes**, so anything downstream that rescales
the fODF must run *after* it, not before. This forced the cap to precede the
now-retired anisotropy modulation, and it applies equally to any peak truncation
you add yourself.

## Where to look next

- Full documentation: root [`README.md`](../../README.md), "fODF outlier capping".
- Measurement report:
  [`Reports/REPORT_fODF_outlier_cap.md`](../../Reports/REPORT_fODF_outlier_cap.md)
- Tests, no data needed:
  [`tests/test_fODF_outlier_cap.m`](../../tests/test_fODF_outlier_cap.m) (the
  method in isolation) and
  [`tests/test_SMI_outlier_cap.m`](../../tests/test_SMI_outlier_cap.m) (the flag
  through `SMI.fit`)

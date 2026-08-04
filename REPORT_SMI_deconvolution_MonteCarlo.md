# Constrained SMI deconvolution vs MSMT-CSD and SSST-CSD

Monte Carlo comparison of spherical deconvolutions on crossing white matter
fibres, in the design of Jeurissen et al. (2014).

**CSD and MSMT-CSD are run by MRtrix3 3.0.4 itself** — `dwi2response`,
`dwi2fod` and `sh2peaks`, the binaries, not a reimplementation of them. The SMI
fODF goes through the same `sh2peaks`, so peak extraction is identical for
every arm and is not something this package implements.

**Everything in this report is simulation.** No number here has touched a
patient scan, or any real scan. Section 2 says exactly which part of the
intended design that compromises and how much it costs.

Code: `deconv_comparison/`. Figures: `fodf_deconv_montecarlo.png`,
`fodf_response_shview.png`. Full result tables: `deconv_tables.md`.

## Findings, shortest form

1. **No method dominates; the three sit on one accuracy-robustness curve.**
   SSST-CSD is the sharpest and the most fragile, MSMT-CSD the most stable and
   the most biased, constrained SMI closest to being both. (section 5)
2. **MSMT-CSD does not resolve a 45 degree crossing in any of 10,000
   realisations, at any SNR** — and that survives every control: exact response
   instead of estimated (0.0% either way), Lmax 8 instead of 6 (0.0% either
   way). Its 60 degree error *is* response-limited: 6.85 degrees with the
   `dwi2response dhollander` response, 2.36 with an exact one. Single-shell CSD
   on the same data resolves 45 degrees 96.9% of the time. (sections 5.3, 7)
3. **Constrained SMI has the best 60 degree accuracy of the three from SNR 50
   down to SNR 10**, and at SNR 10 still returns the correct fibre count in
   99.8% of 60 degree crossings where SSST-CSD manages 49.6% and puts 0.63
   spurious peaks in every voxel. (sections 5.2-5.4)
4. **`lambda_nonneg` is now 1, not 10** (`SMI.m:968`). At 10 the constrained
   deconvolution never resolved a 45 degree crossing either; at 1 it resolves
   81.1% at SNR 50. (section 6)
5. **SMI's SH basis is MRtrix's exactly, but only at `CS_phase = 0`.** At the
   default `CS_phase = 1` the two differ by `(-1)^m`, which is a 180 degree
   rotation about z of every fODF — verified with MRtrix's own `sh2amp` and
   `sh2peaks`. (section 4.1)
6. **Every estimated response is 15-40% blunter than the delta response of the
   kernel that generated the data**, by all three `dwi2response` algorithms.
   Handing CSD an exact response, as earlier work in this repository did, is
   idealised rather than generous. (section 2)

---

## 1. What was asked and what was done

| asked | done |
|---|---|
| estimate realistic response functions from real data with `dwi2response dhollander` | `dwi2response dhollander` is run, the MRtrix binary — but on a synthetic phantom, because there is no real data in this repository |
| also test older approaches | `dwi2response tournier` and `dwi2response fa` estimated alongside |
| noise-free signal vectors for 15, 45 and 60 degree crossings by forward convolution | done, plus a single fibre reference |
| complex Gaussian noise, Rician magnitude, SNR 5 to infinity | SNR 5, 10, 20, 30, 50 and a noise free arm |
| 10,000 realisations | 10,000 per condition per SNR |
| compare constrained SMI to CSD and MSMT-CSD; bias, std, spurious peaks | three arms, all peaks from `sh2peaks` |

**Division of labour.** `gen_montecarlo.m` (Octave/MATLAB) synthesises the
signal, fits SMI, and writes both the DWI and the SMI fODF as MRtrix images.
`run_mrtrix.sh` runs MRtrix. `score_mrtrix.py` reads what MRtrix wrote and does
the bookkeeping. Nothing crosses that line: the only MRtrix behaviour this
package implements is reading and writing its image format (`mrtrix_io.m`,
`mrtrix_io.py`), and that is verified against `mrinfo` and `mrconvert`.

An earlier version of this comparison used dipy for CSD and MSMT-CSD and a
reimplementation of `dwi2response dhollander`. Its conclusions about MSMT-CSD
reproduce here under the real thing; its numbers do not appear in this report.

---

## 2. Where the response functions come from, and what that costs

The responses are estimated by `dwi2response dhollander` (for MSMT-CSD) and
`dwi2response tournier` (for single-shell CSD), the MRtrix commands, with
`-erode 0` because the phantom has no brain edge and `-lmax 0,6,6,6` to hold
every arm at the same angular order.

**But they are estimated from a synthetic phantom, not a brain.**
`gen_phantom.m` builds 9216 voxels from SMI's own forward model: single-fibre
WM with dispersion spread over kappa 10-26, two- and three-fibre crossings,
WM/GM and WM/CSF partial volume, grey matter, CSF, 10% per-voxel jitter on
every kernel parameter, at SNR 30.

What that does buy:

* the responses are **blunted by dispersion** the way real ones are. At
  b = 3 ms/µm², normalised to `l = 0`:

  | | `r_2/r_0` | `r_4/r_0` | `r_6/r_0` |
  |---|---|---|---|
  | `dwi2response dhollander` | -0.703 | 0.342 | -0.107 |
  | `dwi2response tournier` | -0.714 | 0.349 | -0.116 |
  | `dwi2response fa` | -0.718 | 0.352 | -0.121 |
  | an exact delta response from the same kernel | -0.829 | 0.442 | -0.179 |

  Every estimator lands 15-40% short of the delta at every order, and the three
  agree with each other to within 2%. Handing CSD the exact kernel response is
  therefore **idealised, not generous**.
* the responses carry the sampling noise of a finite voxel selection.

What it does **not** buy: any statement about how the selection behaves on real
tissue. The phantom's crossings are at 40-90 degrees, which a single-fibre
criterion rejects easily; on real data the interesting failures are partial
volume, pathology, and crossings near the resolution limit, none of which this
phantom stresses.

**One thing to keep in mind when reading section 5**: because the same forward
model generates both the phantom and the Monte Carlo signals, CSD and MSMT-CSD
face a *milder* model mismatch here than they would on real data, where the
true kernel is not a Standard Model kernel at all. SMI faces none. That is the
largest structural advantage SMI has in this experiment and nothing here
corrects for it.

---

## 3. The simulation

**Protocol.** Three shells at b = 1, 2, 3 ms/µm² with 90 electrostatically
repulsed directions each, plus 18 b = 0 volumes: 288 volumes, the acquisition
MSMT-CSD was introduced on. Written into the MRtrix header as a `dw_scheme` in
s/mm², so MRtrix reads the gradient table from the image itself.

**Ground truth.** Two equal Watson populations, kappa = 16, at 0 (single
fibre), 15, 45 and 60 degrees, convolved with a white matter kernel
`[f Da Depar Deperp fw] = [0.60 2.0 2.0 0.50 0.02]` through SMI's forward
model. The truth is expanded to **Lmax 8**, not to the Lmax 6 every method fits
at, so no method is handed a band-limited target it can represent exactly. The
truth is also written as an SH image and run through `sh2peaks`, which is where
the "band limited truth" column in section 5.2 comes from.

**Noise.** `sqrt((S + sigma*n1)^2 + (sigma*n2)^2)` with `n1, n2` standard
normal and `sigma = 1/SNR`: complex Gaussian noise, magnitude taken, i.e.
Rician. S0 = 1 by construction. 10,000 independent realisations per condition
per SNR; the noise-free arm uses `SNR = 1e4` and 500 realisations.

**The three deconvolutions**, all at Lmax 6, all on the same image:

| arm | what |
|---|---|
| SMI constrained | `SMI.fit` with `fODF_regularization.flag_nonneg = 1` at the shipped defaults (`lambda_nonneg = 1`, `lambda_tikhonov = 0.3`), `CS_phase = 0` |
| SSST-CSD | `dwi2fod csd` on the b = 3 shell, `dwi2response tournier` response |
| MSMT-CSD | `dwi2fod msmt_csd` on all shells, `dwi2response dhollander` 3-tissue responses |

SMI estimates a **kernel per voxel** from the data itself; CSD and MSMT-CSD get
**one global response**. That asymmetry is not a flaw in the comparison, it is
the methodological difference being compared.

**Scoring.** `sh2peaks -num 4` on every arm, then one post hoc rule applied
identically: a peak counts if its *anisotropic* amplitude (amplitude minus the
fODF's isotropic term `c_00 Y_00`) is positive and at least 0.30 of the largest
anisotropic amplitude in that voxel. Subtracting the isotropic term is what
makes the threshold mean the same thing for SMI, whose fODF has unit mass and a
fixed `1/(4*pi)` floor, and for CSD, whose `l = 0` term is its apparent fibre
density.

---

## 4. What was verified before anything was measured

Comparing four fODFs computed by three different codebases only means something
if the conventions line up. `deconv_comparison/check_conventions.py` checks
them, and `test_SMI_response_helpers.m` checks the MATLAB side (8 tests, all
passing).

| check | result |
|---|---|
| SMI's SH basis vs dipy's `descoteaux07`: is there an exact linear map? | yes, residual `3.9e-15`. It is a **signed permutation**, not the identity |
| a delta along a known axis, mapped and run through the peak finder | recovered at that axis, `0.0000` deg |
| the convolution rule `s_lm = c_lm r_l sqrt(4 pi/(2l+1))` vs SMI's forward model | `5.6e-16` |
| the MRtrix response file written by `example_SMI_response_shview.m` vs `kernel.py` | `3.8e-9`, i.e. the `%.8g` write precision |
| the zonal response `r_l = K_l sqrt((2l+1) 4 pi)` vs SMI's forward model, in MATLAB | `1.2e-15` |

### 4.1 SMI's basis is MRtrix's at `CS_phase = 0`, and rotated at `CS_phase = 1`

This was an open item in `README for Claude` section 8 item 4 — "verify SMI's SH
basis against MRtrix's, every peak orientation shipped to `tckgen` depends on
it". It is now measured **against MRtrix itself**, not against a
reimplementation of it: `deconv_comparison/check_mrtrix_basis.sh` hands
`sh2amp` an SH image whose 28 voxels each hold one unit basis vector, so
`sh2amp`'s output *is* MRtrix's basis matrix, and compares it against
`SMI.get_even_SH`.

| | map from MRtrix's basis to SMI's | residual |
|---|---|---|
| `options.CS_phase = 0` | **the identity** | `3.3e-8` (the float32 image round trip) |
| `options.CS_phase = 1` (**the default**) | `diag((-1)^m)` | `3.3e-8` |

Off-diagonal terms are exactly zero in both cases.

**`CS_phase = 0` means SMI's `out.plm` needs no conversion at all** — the
coefficients are already MRtrix's. That is what this package now uses
throughout, so nothing sits between the SMI arm and the MRtrix arms.

`CS_phase = 1` is a different story. For even `l`, multiplying real SH
coefficients by `(-1)^m` is exactly a **180 degree rotation about the z axis**.
Measured, with MRtrix's own `sh2peaks`, on a sharp fODF along
`[0.3010 -0.5017 0.8127]`:

| written with | `sh2peaks` returns | error |
|---|---|---|
| `CS_phase = 0` | `[ 0.3006 -0.5010  0.8116]` | **0.00 deg** |
| `CS_phase = 1` | `[-0.3006  0.5010  0.8116]` | **71.50 deg** — `x` and `y` negated |

So: **an SMI fODF written out with the default `CS_phase = 1` and tracked in
MRtrix is a fibre field rotated 180 degrees about z.** Set
`options.CS_phase = 0` before writing `out.plm` for MRtrix, or apply `(-1)^m`
on the way out.

> **Correction.** An earlier version of this section reported the same
> discrepancy as a *reflection through `z = 0`*, and reported an extra
> `sqrt(2)` normalisation at `m != 0`. Both were artefacts of measuring
> against dipy's `tournier07` basis instead of MRtrix's — dipy's is not
> MRtrix's. The magnitude of the error (~72 degrees) was right; the geometry
> and the normalisation were not.

---

## 5. Results

Figure: `fodf_deconv_montecarlo.png`. Every table below is a slice of
`deconv_tables.md`, generated by `deconv_comparison/tables.py` from the peaks
MRtrix's `sh2peaks` returned — no number in this report is transcribed by hand.

### 5.1 In one paragraph

The three methods are three points on the same accuracy-robustness curve, and
none of them dominates. **SSST-CSD is the sharpest** — at SNR 50 it resolves 45
degree crossings in 96.9% of realisations at 2.84 degrees, better than anything
else here — **and the first to fall apart**: by SNR 5 it returns the right
number of fibres in 15.6% of single-fibre voxels and puts 1.6 spurious peaks in
every one. **MSMT-CSD is the most stable and the most biased**: essentially no
spurious peaks at any SNR, an angular error that barely moves, and a 45 degree
crossing it never resolves in 10,000 realisations at any SNR. **Constrained SMI
sits between them and is closest to being both**: within 1.4 degrees of
SSST-CSD's 45 degree accuracy at high SNR, the best 60 degree accuracy of the
three from SNR 50 down to 10, and at SNR 10 it still resolves 99.8% of 60 degree
crossings where SSST-CSD manages 49.6%.

### 5.2 Angular error of the largest peak, degrees (median)

_single fibre_ (band limited truth: 0.00°)

| SNR | SMI constrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|
| noise free | **0.03** | 0.07 | 0.32 |
| 50 | **0.27** | 0.46 | 0.35 |
| 30 | 0.45 | 0.74 | **0.42** |
| 20 | 0.67 | 1.10 | **0.57** |
| 10 | 1.38 | 2.24 | **1.08** |
| 5 | 3.03 | 5.12 | **2.33** |

_crossing 45°_ (band limited truth: 1.64°)

| SNR | SMI constrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|
| noise free | 6.86 | **2.47** | 22.34 |
| 50 | 6.52 | **2.84** | 22.25 |
| 30 | 6.47 | **3.50** | 22.07 |
| 20 | 6.69 | **4.71** | 21.82 |
| 10 | 9.00 | **8.35** | 20.19 |
| 5 | **12.28** | 13.40 | 16.37 |

_crossing 60°_ (band limited truth: 1.30°)

| SNR | SMI constrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|
| noise free | **1.50** | 1.62 | 7.21 |
| 50 | **1.43** | 1.61 | 6.85 |
| 30 | **1.45** | 1.87 | 6.50 |
| 20 | **1.70** | 2.49 | 6.22 |
| 10 | **3.25** | 5.08 | 5.29 |
| 5 | 7.97 | 11.92 | **7.60** |

MSMT-CSD's ~22 degrees at 45 degrees is not noise. It is what the largest peak
of a single unresolved lobe measures against either true fibre, and it is
within 0.5 degrees of its own noise-free value at every SNR down to 20 — the
failure is structural, and section 7 shows it survives both an exact response
and Lmax 8.

The 15 degree crossing (`deconv_tables.md`) is unresolvable at Lmax 6 by
construction: the band limited ground truth itself presents a single peak 7.50
degrees from each fibre, and all three methods reproduce that within 0.3
degrees at SNR 50. What that column measures is how gracefully each degrades.

### 5.3 Realisations returning the correct number of fibres, %

_crossing 45°_

| SNR | SMI constrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|
| noise free | 100.0 | 100.0 | 0.0 |
| 50 | 81.1 | **96.9** | 0.0 |
| 30 | 63.4 | **83.5** | 0.0 |
| 20 | 51.5 | **68.3** | 0.0 |
| 10 | 36.5 | **40.7** | 0.0 |
| 5 | **43.9** | 4.2 | 12.2 |

_crossing 60°_

| SNR | SMI constrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|
| noise free | 100.0 | 100.0 | 100.0 |
| 30 | **100.0** | 100.0 | 99.6 |
| 20 | **100.0** | 100.0 | 96.8 |
| 10 | **99.8** | 49.6 | 87.2 |
| 5 | 45.5 | 1.6 | **82.7** |

_single fibre_ — all three are at 100.0% down to SNR 10. At SNR 5 SMI is at
96.8% and MSMT-CSD at 100.0%, while **SSST-CSD collapses to 15.6%**: five in
six single-fibre voxels come back with more than one peak.

The 45 degree column at SNR 5 is a trap worth naming. SMI's 43.9% and
MSMT-CSD's 12.2% are not recoveries — they are noise splitting a single lobe
into two, which happens to be the right *count*. The angular errors in 5.2
(12.28 and 16.37 degrees) say the directions are wrong.

### 5.4 Spurious peaks per realisation

Mean number of peaks above the true count.

| condition | SNR | SMI constrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|---|
| single fibre | 10 | **0.000** | 0.007 | **0.000** |
| single fibre | 5 | 0.033 | 1.602 | **0.000** |
| crossing 45° | 10 | **0.001** | 0.452 | **0.000** |
| crossing 45° | 5 | 0.316 | 1.734 | **0.001** |
| crossing 60° | 10 | 0.002 | 0.633 | **0.001** |
| crossing 60° | 5 | 0.670 | 1.845 | **0.016** |

Down to SNR 10, constrained SMI and MSMT-CSD are both at or below 0.002 per
voxel while SSST-CSD is at 0.45-0.63 — that is a spurious direction in roughly
every second voxel, for a tractography algorithm to follow. At SNR 5 SMI's
count rises to 0.3-0.7 and MSMT-CSD's stays at zero.

### 5.5 Angular correlation coefficient against the band limited truth

Scale free and uses the whole fODF rather than its maxima, so it is the one
metric on which all three are directly comparable.

_crossing 60°_ (mean over 10,000 realisations)

| SNR | SMI constrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|
| noise free | 0.9753 | **0.9903** | 0.8448 |
| 50 | 0.9736 | **0.9848** | 0.8470 |
| 30 | 0.9708 | **0.9737** | 0.8487 |
| 20 | **0.9646** | 0.9488 | 0.8497 |
| 10 | **0.9190** | 0.8100 | 0.8518 |
| 5 | 0.7213 | 0.4805 | **0.7846** |

The crossover is at about SNR 20-30. Above it SSST-CSD reconstructs the fODF
best; below it constrained SMI does, and by SNR 5 MSMT-CSD's flat curve has
overtaken both. The same ordering and the same crossover appear at 45 degrees,
at 15 degrees and on the single fibre (`deconv_tables.md`).

### 5.6 Peak amplitude

Normalised to each method's own noise-free median, so methods whose fODFs are
not on the same scale can still be compared. Full table in `deconv_tables.md`.
The pattern is the one that matters for any amplitude-thresholded tractography:
**noise inflates the peak of the sharper methods and leaves MSMT-CSD's almost
alone.** At SNR 5 on a 45 degree crossing SSST-CSD's peak is well above its own
noise-free value with a coefficient of variation near 0.2, while constrained
SMI moves less and MSMT-CSD least. A fixed `-cutoff` therefore means a
different thing at different SNR, in a method-dependent direction.

---

## 6. Why `lambda_nonneg` is now 1

`README for Claude` section 2 item 3 recorded that the shipped
`lambda_nonneg = 10` was what closed the 45 degree crossing, and asked for a
cross-check against an independent metric before changing anything. This is
that cross-check: `deconv_comparison/sweep_nonneg.m`, 2000 realisations per
condition at SNR 30, one `SMI.fit` and then one deconvolution per setting from
the same kernel, scored by the same `sh2peaks` as the rest of the report.

**Realisations returning the correct number of fibres, %**

| setting | single | 15° | 45° | 60° |
|---|---|---|---|---|
| non-negativity off | 100.0 | 0.0 | **95.9** | 99.8 |
| `lambda_nonneg = 1` (**now the default**) | 100.0 | 0.0 | 55.4 | 100.0 |
| `lambda_nonneg = 3` | 100.0 | 0.0 | 0.2 | 100.0 |
| `lambda_nonneg = 10` (the old default) | 100.0 | 0.0 | 0.0 | 100.0 |

**Angular error of the largest peak, degrees (median)**

| setting | single | 15° | 45° | 60° |
|---|---|---|---|---|
| non-negativity off | 1.11 | 6.78 | **2.81** | 2.40 |
| `lambda_nonneg = 1` | 0.43 | 7.24 | 7.27 | 1.37 |
| `lambda_nonneg = 3` | **0.34** | 7.29 | 18.30 | **0.87** |
| `lambda_nonneg = 10` | 0.35 | 7.30 | 21.28 | 1.94 |

**Angular correlation coefficient against the band limited truth (mean)**

| setting | single | 15° | 45° | 60° |
|---|---|---|---|---|
| non-negativity off | 0.9751 | 0.9707 | 0.9425 | 0.9448 |
| `lambda_nonneg = 1` | **0.9800** | **0.9857** | **0.9656** | **0.9718** |
| `lambda_nonneg = 3` | 0.9572 | 0.9689 | 0.9311 | 0.9439 |
| `lambda_nonneg = 10` | 0.9298 | 0.9465 | 0.8869 | 0.9032 |

**Spurious peaks per realisation**: exactly 0.000 for every constrained
setting at every condition. With the constraint **off** they appear — 0.026 at
45 degrees and 0.003 at 60 degrees. That is what non-negativity buys, and the
sweep says it is bought in full at `lambda_nonneg = 1`; everything above 1 pays
more for nothing.

Three conclusions, and the first is why the default changed.

1. **`lambda_nonneg = 1` has the best angular correlation at every condition**
   — including against turning the constraint off entirely — and it already
   buys all of the spurious-peak suppression. Nothing above it improves any
   metric except the 60 degree median error at 3 (0.87 vs 1.37 degrees), which
   costs 55.2 percentage points of 45 degree resolution. **The default is now
   1** (`SMI.m:968`).
2. **The old default of 10 was on nobody's Pareto front.** It is worse than 3
   on 45 degrees, on 60 degrees and on every ACC, and worse than 1 on
   everything.
3. **Tikhonov still does nothing.** `lambda_tikhonov` 0.3 versus 0 at
   `lambda_nonneg = 10` moves the 45 degree error from 21.28 to 21.27 degrees
   and the ACC in the fourth decimal. That confirms `README for Claude`
   section 2 item 1 on a third independent metric.

One honest caveat about the direction of this change. At `lambda_nonneg = 1`
the 45 degree crossing is resolved 55.4% of the time here and 63.4% in the
10,000-realisation arm of section 5 — better than 10, far short of turning the
constraint off (95.9%). The reason to stop at 1 rather than 0 is the spurious
peak column and the SNR 5 behaviour in section 5.4, not the 45 degree column.

---

## 7. Controls: is MSMT-CSD being treated fairly?

A comparison in which MSMT-CSD never resolves a 45 degree crossing is equally
consistent with "MSMT-CSD was set up badly here", so it was run twice more at
SNR 50, changing one thing at a time.

**Control 1: the exact response.** Both CSD variants rerun with the delta-fODF
response of the very kernel the signals were convolved with — better than any
estimator could produce (`write_exact_response.py`, `run_mrtrix.sh control`).

**Control 2: Lmax 8.** Both rerun at `-lmax 8`, which is what MRtrix would
choose by itself on 288 volumes; the headline comparison holds every arm at
Lmax 6 so SMI is not disadvantaged.

| | 45° resolved | 45° error | 60° error | single fibre error |
|---|---|---|---|---|
| SSST-CSD, `tournier` response, Lmax 6 | 96.9% | 2.84° | 1.61° | 0.46° |
| SSST-CSD, **exact** response, Lmax 6 | **99.4%** | **1.73°** | 1.59° | 0.54° |
| SSST-CSD, `tournier` response, **Lmax 8** | 99.9% | **1.14°** | **1.11°** | 0.44° |
| MSMT-CSD, `dhollander` responses, Lmax 6 | **0.0%** | 22.25° | 6.85° | 0.35° |
| MSMT-CSD, **exact** responses, Lmax 6 | **0.0%** | 22.04° | **2.36°** | 0.22° |
| MSMT-CSD, `dhollander` responses, **Lmax 8** | **0.0%** | 21.10° | **2.50°** | **0.25°** |

Three things follow.

* **MSMT-CSD's 45 degree failure is neither a response problem nor a band-limit
  problem.** 0.0% of 10,000 realisations in all three configurations. The band
  limited ground truth resolves it 100% of the time at Lmax 6, so the
  information is there and MSMT-CSD's fit does not use it.
* **MSMT-CSD's 60 degree error is a response problem, and also a band-limit
  problem.** Either fix takes it from 6.85 degrees to about 2.4. Both effects
  are real and roughly the same size.
* **Single-shell CSD's 45 degree accuracy is response-limited**, and it gets
  most of the way there on the estimated response alone: 2.84 degrees with
  `dwi2response tournier` against 1.73 with the exact one, which is within 0.1
  degrees of the band limited ground truth's own 1.64.

The mechanism for the 45 degree failure is not measured here, but the shape of
the result points at one: MSMT-CSD fits all three shells jointly, and b = 1 and
b = 2 carry very little angular contrast at `l >= 4` for a 45 degree crossing,
so they pull the joint solution towards a single lobe. Single-shell CSD sees
only b = 3, which is where the contrast is. If that is right, the practical
reading is that **MSMT-CSD's multi-shell advantage is tissue separation, not
angular resolution** — and this experiment deliberately contains no tissue
heterogeneity for it to separate (section 8).

---

## 8. What this does not show

Worth being explicit about, because several of these are easy to read into the
tables and none of them is measured here.

* **Nothing has touched real data**, including the "real data" the responses
  come from (section 2). The single most valuable next step is to run
  `dwi2response dhollander` on an actual brain and compare the response it
  returns against the phantom-derived ones here.
* **Only one microstructure.** Every Monte Carlo voxel has
  `[f Da Depar Deperp fw] = [0.60 2.0 2.0 0.50 0.02]`. SMI estimates a kernel
  per voxel and CSD does not; a population with varying microstructure is where
  that difference should pay off, and it is not tested here.
* **Only two equal fibre populations.** No unequal volume fractions, no three-way
  crossings, no fanning. Jeurissen et al. (2014) swept the volume fraction ratio
  as well as the angle; this does not.
* **No partial volume with grey matter or CSF in the Monte Carlo conditions.**
  That is precisely where MSMT-CSD is supposed to win, and it is deliberately
  outside this comparison — `REPORT_fODF_freewater.md` covers it. Reading these
  tables as "MSMT-CSD is worse than SSST-CSD" would be exactly the wrong
  conclusion: this experiment contains none of the tissue heterogeneity MSMT-CSD
  exists to handle.
* **One protocol, one Lmax for the headline.** Three shells, 90 directions
  each, Lmax 6 everywhere; section 7 shows what Lmax 8 changes for the CSD
  arms, and it is not nothing. The 15 degree crossing is unresolvable at Lmax 6 by construction —
  the band limited ground truth itself presents a single peak — so that column
  measures how gracefully each method degrades, not whether it can resolve 15
  degrees.
* **The comparison is not blinded to the forward model.** Both the phantom and
  the Monte Carlo signals come from SMI's own Standard Model forward model, so
  CSD and MSMT-CSD face a milder model mismatch than they would on real tissue,
  and SMI faces none at all. This is the largest structural advantage SMI has in
  this experiment and no number here corrects for it.

---

## 9. Reproducing it

```
cd deconv_comparison
python3 setup_protocol.py                                  # protocol + eval sphere
octave --eval "oct_path; gen_phantom(30,'p30')"            # response phantom
./run_mrtrix.sh responses                                  # dwi2response x3

for s in 50 30 20 10 5; do
  octave --eval "oct_path; gen_montecarlo($s,10000,'snr$s')"   # signal + SMI
  ./run_mrtrix.sh fit snr$s                                    # dwi2fod + sh2peaks
done
octave --eval "oct_path; gen_montecarlo(1e4,500,'nf')"
./run_mrtrix.sh fit nf

./check_mrtrix_basis.sh                                    # SH basis vs MRtrix
python3 tables.py nf snr50:50 snr30:30 snr20:20 snr10:10 snr5:5 \
    | tee ../deconv_tables.md
python3 figure_mc.py ../fodf_deconv_montecarlo.png \
    snr50:50 snr30:30 snr20:20 snr10:10 snr5:5
python3 figure_response.py p30 ../fodf_response_shview.png

octave --eval "oct_path; sweep_nonneg(30,2000,'sw30')"     # section 6
./run_mrtrix.sh sweep sw30 6
python3 score_sweep.py sw30 30
./run_mrtrix.sh control snr50                              # section 7
```

`run_all.sh` does all of it. The SNR arms are independent and were run three or
four at a time on four cores; one arm is about 12 minutes of `SMI.fit` and
under a minute of MRtrix. The whole thing is about an hour of wall clock on
four cores, plus about 25 minutes for the sweep and controls.

Environment: MRtrix3 **3.0.4**, GNU Octave 8.4 with `statistics` and three
shims for MATLAB functions Octave lacks (`deconv_comparison/stubs/`), Python
3.11 with `numpy` and `matplotlib` for the bookkeeping and figures.

---

## 10. What to do next, in order of value

1. **Run `dwi2response dhollander` on a real brain** and compare the responses
   against the phantom-derived ones in section 2. That is the one substitution
   left in this package, and removing it is cheap.
2. **Find out why MSMT-CSD never resolves 45 degrees** (section 7). The
   shell-weighting explanation is a hypothesis, not a measurement; refitting
   MSMT-CSD on b = 0 and b = 3 alone would test it in one command. If it holds,
   it is a statement about multi-shell deconvolution worth making carefully,
   and if it does not, something else is going on and this report should say so.
3. **Add unequal fibre volume fractions and three-way crossings.** The
   machinery takes them — `gen_montecarlo.m` builds its conditions from a list
   of axes — and they are the part of the Jeurissen design that is missing.
4. **Re-run with a heterogeneous microstructure per voxel.** SMI's per-voxel
   kernel should pay off there and cannot pay off in the current design, where
   every voxel has the same kernel.
5. **Put grey matter and CSF partial volume back in.** This experiment removes
   exactly the thing MSMT-CSD exists for, which makes its numbers here easy to
   misread. `REPORT_fODF_freewater.md` covers part of it; a combined design
   would be better.

---

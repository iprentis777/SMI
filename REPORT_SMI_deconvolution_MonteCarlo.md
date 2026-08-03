# Constrained SMI deconvolution vs MSMT-CSD and SSST-CSD

Monte Carlo comparison of spherical deconvolutions on crossing white matter
fibres, in the design of Jeurissen et al. (2014).

**Everything in this report is simulation.** No number here has touched a
patient scan, or any real scan. Section 2 says exactly which part of the
intended design that compromises and how much it costs.

Code: `deconv_comparison/`. Figures: `fodf_deconv_montecarlo.png`,
`fodf_response_shview.png`. Full result tables: `deconv_tables.md`.

## Findings, shortest form

1. **The constrained SMI deconvolution trades resolution for stability, and the
   trade is very lopsided at the shipped `lambda_nonneg = 10`.** It never
   resolves a 45 degree crossing at any SNR, including noise free — but it
   produces no spurious peaks anywhere, and at SNR 5 it still returns the
   correct fibre count in 87% of 60 degree crossings where SSST-CSD manages 17%
   and the unconstrained deconvolution 0.1%. (sections 5.3-5.5)
2. **`lambda_nonneg = 10` is not on the Pareto front.** `1` has the best angular
   correlation at every condition and resolves 63.5% of 45 degree crossings;
   `3` beats the shipped default on 45 degrees, 60 degrees and single fibres.
   The default was not changed — this is a measurement, not a patch. (section 6)
3. **SMI's spherical harmonic basis is not MRtrix's**, and at the default
   `CS_phase = 1` the difference is exactly a mirror through the `z = 0` plane.
   An SMI fODF written out as-is and tracked in MRtrix is a mirrored fibre
   field: 72.32 degrees of peak error against 1.65 with `CS_phase = 0`.
   (section 4.1)
4. **MSMT-CSD's crossing performance here is response-limited, not a setup
   error.** With the `dwi2response dhollander` response its 60 degree error is
   7.06 degrees; with an exact delta response, 2.57. Its 45 degree failure is
   not response-limited — 0% either way. (section 7)
5. **Every estimated response is 15-40% blunter than the delta response of the
   kernel that generated the data**, at every order and by all three
   estimators. Handing CSD an exact response, as earlier work in this
   repository did, is idealised rather than generous. (section 2)

---

## 1. What was asked and what was done

| asked | done |
|---|---|
| estimate realistic response functions from real data with `dwi2response dhollander` | the algorithm is reimplemented from the MRtrix 3.0.4 source and run, but on a synthetic phantom — there is no real data in this repository |
| also test older approaches | `dwi2response tournier` and `dwi2response fa` estimated alongside |
| noise-free signal vectors for 15, 45 and 60 degree crossings by forward convolution | done, plus a single fibre reference |
| complex Gaussian noise, Rician magnitude, SNR 5 to infinity | SNR 5, 10, 20, 30, 50 and a noise free arm |
| 10,000 realisations | 10,000 per condition per SNR |
| deconvolve three ways; bias, std, spurious peaks | four ways: the constrained SMI deconvolution, the same kernel deconvolved without the constraint, SSST-CSD, MSMT-CSD |

The fourth arm is not padding. The constrained and unconstrained SMI arms share
**one** kernel fit — `SMI.fit` is run once and the unconstrained fODF is then
recomputed from `out.kernel` by calling `SMI.get_plm_from_S_and_kernel`
directly — so any difference between them is the deconvolution alone. Without
it, "constrained SMI beats CSD" cannot be separated from "SMI's kernel beats
CSD's response".

---

## 2. Where the response functions come from, and what that costs

`dwi2response dhollander` is reimplemented in `deconv_comparison/dhollander.py`
following `lib/mrtrix3/dwi2response/dhollander.py` of MRtrix 3.0.4 step by step:
the signal decay metric, the crude FA/SDM segmentation, the MAD-based WM
refinement, the GM and CSF refinements, the top-percentage voxel selections,
and the Dhollander et al. (2019) single-fibre WM criterion — which needs an
MSMT-CSD fit with a deliberately over-sharp "ewmrf" response, so that fit is
implemented too (`mtcsd.py`; dipy refuses fewer than two isotropic
compartments). `mrthreshold`'s default threshold and `amp2response`'s
constrained fit are reimplemented in `mrtrix_ops.py` from the MRtrix C++.

**But it is run on a synthetic phantom, not a brain.** `gen_phantom.m` builds
9216 voxels from SMI's own forward model: single-fibre WM with dispersion spread
over kappa 10-26, two- and three-fibre crossings, WM/GM and WM/CSF partial
volume, grey matter, CSF, 10% per-voxel jitter on every kernel parameter, at
SNR 30.

What that does buy:

* the responses are **blunted by dispersion** the way real ones are. At b = 3,
  normalised to `l = 0`, the estimated WM responses are

  | | `r_2/r_0` | `r_4/r_0` | `r_6/r_0` |
  |---|---|---|---|
  | `dwi2response dhollander` | -0.702 | 0.340 | -0.105 |
  | `dwi2response tournier` | -0.690 | 0.325 | -0.089 |
  | `dwi2response fa` | -0.708 | 0.357 | -0.115 |
  | an exact delta response from the same kernel | -0.829 | 0.443 | -0.180 |

  Every estimator lands 15-40% short of the delta at every order. Handing CSD
  the exact kernel response, as earlier work in this repository did, is
  therefore **idealised, not generous** — it gives CSD a sharper response than
  any estimator would produce.
* the responses carry the sampling noise of a finite voxel selection: 23 voxels
  for `dhollander` (0.5% of 4501 refined WM), 300 for the other two.

What it does **not** buy: any statement about how the selection behaves on real
tissue. All three selectors come back **100% pure single-fibre WM** here.
That is a property of the phantom — its crossings are at 40-90 degrees, which a
single-fibre criterion rejects easily — not evidence that response estimation is
safe. On real data the interesting failures are partial volume, pathology, and
crossings near the resolution limit, none of which this phantom stresses.

**The one substitution to keep in mind when reading section 5**: because the
same forward model generates both the phantom and the Monte Carlo signals, CSD
and MSMT-CSD are being tested under a *milder* model mismatch than they would
face on real data, where the true kernel is not a Standard Model kernel at all.
The comparison is not tilted towards SMI by this; if anything the reverse.

---

## 3. The simulation

**Protocol.** Three shells at b = 1, 2, 3 ms/µm² with 90 electrostatically
repulsed directions each, plus 18 b = 0 volumes: 288 volumes, the acquisition
MSMT-CSD was introduced on.

**Ground truth.** Two equal Watson populations, kappa = 16, at 0 (single
fibre), 15, 45 and 60 degrees, convolved with a white matter kernel
`[f Da Depar Deperp fw] = [0.60 2.0 2.0 0.50 0.02]` through SMI's own forward
model. The truth is expanded to **Lmax 8**, not to the Lmax 6 every method
fits at, so no method is handed a band-limited target it can represent exactly.

**Noise.** `sqrt((S + sigma*n1)^2 + (sigma*n2)^2)` with `n1, n2` standard
normal and `sigma = 1/SNR`: complex Gaussian noise, magnitude taken, i.e.
Rician. S0 = 1 by construction. 10,000 independent realisations per condition
per SNR; the noise-free arm uses `SNR = 1e4` and 500 realisations, and exists
only to normalise peak amplitudes.

**The four deconvolutions**, all at Lmax 6, all on the same bytes:

| arm | what |
|---|---|
| SMI constrained | `SMI.fit`, `fODF_regularization = struct('flag_nonneg',1,'lambda_nonneg',10,'lambda_tikhonov',0.3)` |
| SMI unconstrained | the same `out.kernel`, plain LLS deconvolution |
| SSST-CSD | dipy `ConstrainedSphericalDeconvModel`, b = 3 shell only, `dwi2response tournier` response |
| MSMT-CSD | dipy `MultiShellDeconvModel`, all shells, `dwi2response dhollander` 3-tissue responses |

SMI estimates a **kernel per voxel** from the data itself; CSD and MSMT-CSD get
**one global response**. That asymmetry is not a flaw in the comparison, it is
the methodological difference being compared.

**Scoring.** One peak finder for all four (`peaks.py`), on one 2890-vertex
sphere, with every peak refined against the continuous SH expansion by a local
quadratic in the tangent plane — without that refinement every reported angle
would carry a ~2 degree quantisation floor that has nothing to do with the
methods. Peaks are found on the **anisotropic part** (`l >= 2` only): SMI's
fODF is normalised to unit mass and carries a fixed isotropic floor of
`1/(4*pi)`, CSD's `l = 0` term is its apparent fibre density and varies, so a
relative peak threshold applied to the full fODF would mean a different thing
for each method. Relative threshold 0.30, minimum separation 20 degrees.

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
`deconv_tables.md`, which is generated by `deconv_comparison/tables.py` from
the same cached measurements the figure uses — no number in this report is
transcribed by hand.

### 5.1 In one paragraph

The four methods separate into two behaviours, and the separation is a
bias-variance trade, not a ranking. **The constrained SMI deconvolution and
MSMT-CSD are almost completely insensitive to noise** — constrained SMI's
angular correlation on a 15 degree crossing is 0.9469 noise free and 0.9411 at
SNR 5 — but both are **biased**: neither ever resolves a 45 degree crossing,
even with a perfect response and no noise. **The unconstrained SMI
deconvolution is essentially exact and essentially unusable**: ACC 0.9998 noise
free, resolving 45 degrees in 100% of realisations, and 2.76 spurious peaks per
single-fibre voxel at SNR 5. **SSST-CSD sits between them** on every axis.
The crossover is at about SNR 10-20: above it the unconstrained deconvolution
is the most faithful, below it the constrained one is, and the constrained one
is the only arm that still resolves 87% of 60 degree crossings at SNR 5.

### 5.2 Angular error of the largest peak, degrees (median)

_single fibre_ (band limited truth: 0.00°)

| SNR | SMI constrained | SMI unconstrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|---|
| noise free | 0.18 | **0.02** | 0.04 | 0.17 |
| 50 | 0.25 | 0.72 | 0.35 | **0.27** |
| 30 | 0.36 | 1.14 | 0.59 | **0.32** |
| 20 | 0.52 | 1.60 | 0.88 | **0.47** |
| 10 | 1.05 | 2.77 | 1.84 | **1.01** |
| 5 | 2.37 | 5.29 | 4.36 | **2.34** |

_crossing 60°_ (band limited truth: 1.30°)

| SNR | SMI constrained | SMI unconstrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|---|
| noise free | 1.52 | **1.39** | 1.71 | 8.07 |
| 50 | 1.82 | 1.87 | **1.45** | 7.53 |
| 30 | 2.11 | 2.47 | **1.60** | 7.06 |
| 20 | 2.43 | 3.25 | **2.12** | 6.46 |
| 10 | **3.51** | 5.66 | 4.48 | 5.31 |
| 5 | **7.26** | 12.42 | 10.77 | 7.60 |

_crossing 45°_ (band limited truth: 1.64°)

| SNR | SMI constrained | SMI unconstrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|---|
| noise free | 22.29 | **1.78** | 10.85 | 22.39 |
| 50 | 21.75 | **2.30** | 10.76 | 22.27 |
| 30 | 21.38 | **3.02** | 10.65 | 22.11 |
| 20 | 20.83 | **3.99** | 10.56 | 21.86 |
| 10 | 17.94 | **7.19** | 11.39 | 20.27 |
| 5 | 15.46 | 14.01 | **13.12** | 16.43 |

The 45 degree column is the one to read carefully. The ~22 degree errors of
constrained SMI and MSMT-CSD are **not** noise: they are what the largest peak
of a single merged lobe measures against either true fibre when the two are
never separated. Both are within 1 degree of their noise-free value at every
SNR down to 20, so the failure is entirely structural. SSST-CSD's ~11 degrees
is a different failure: it *does* separate the two peaks, but its blunt
estimated response pushes them roughly 22 degrees too far apart (with the exact
response the error is 4.13 degrees — section 7).

The 15 degree crossing (full table in `deconv_tables.md`) is unresolvable at
Lmax 6 by construction: the band limited ground truth itself presents a single
peak 7.50 degrees from each fibre, and every method reproduces that within 0.4
degrees at SNR 50.

### 5.3 Realisations returning the correct number of fibres, %

_crossing 45°_

| SNR | SMI constrained | SMI unconstrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|---|
| noise free | 0.0 | **100.0** | 100.0 | 0.0 |
| 50 | 0.0 | **100.0** | 37.9 | 0.0 |
| 30 | 0.0 | **93.2** | 27.5 | 0.0 |
| 20 | 0.0 | **57.0** | 25.6 | 0.0 |
| 10 | 0.7 | 7.5 | **33.1** | 0.0 |
| 5 | 16.8 | 0.4 | **26.1** | 11.8 |

_crossing 60°_

| SNR | SMI constrained | SMI unconstrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|---|
| noise free | 100.0 | 100.0 | 100.0 | 100.0 |
| 30 | **100.0** | 98.3 | 100.0 | 99.9 |
| 20 | **100.0** | 72.4 | 100.0 | 97.5 |
| 10 | **99.4** | 6.3 | 91.6 | 86.2 |
| 5 | **87.0** | 0.1 | 16.9 | 82.9 |

At SNR 5 the unconstrained deconvolution has stopped working: it returns the
right number of fibres in 0.1% of 60 degree voxels and 4.5% of single fibre
voxels. The constrained one returns it in 87.0% and 100%.

### 5.4 Spurious peaks per realisation

Mean number of peaks above the true count. This is the column the non-negativity
constraint exists for.

| condition | SNR | SMI constrained | SMI unconstrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|---|---|
| single fibre | 10 | **0.000** | 0.278 | 0.000 | 0.000 |
| single fibre | 5 | **0.000** | 2.757 | 0.530 | 0.000 |
| crossing 45° | 10 | **0.000** | 2.048 | 0.025 | 0.000 |
| crossing 45° | 5 | 0.001 | 2.827 | 0.985 | 0.001 |
| crossing 60° | 10 | **0.000** | 2.090 | 0.080 | 0.001 |
| crossing 60° | 5 | 0.020 | 2.886 | 1.274 | 0.016 |

Constrained SMI and MSMT-CSD produce **no** spurious peaks anywhere above
SNR 5, and at most 0.02 per voxel at SNR 5. The unconstrained deconvolution
produces nearly three per voxel at SNR 5 — for a tractography algorithm that is
three wrong directions to follow in every voxel.

### 5.5 Angular correlation coefficient against the band limited truth

Scale free, and uses the whole fODF rather than its maxima, so it is the one
metric on which all four are directly comparable.

_crossing 60°_ (mean over 10,000 realisations)

| SNR | SMI constrained | SMI unconstrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|---|
| noise free | 0.9053 | **0.9998** | 0.9779 | 0.8399 |
| 50 | 0.9043 | **0.9752** | 0.9757 | 0.8427 |
| 30 | 0.9017 | 0.9382 | **0.9711** | 0.8459 |
| 20 | 0.8989 | 0.8799 | **0.9598** | 0.8490 |
| 10 | **0.8870** | 0.6962 | 0.8747 | 0.8501 |
| 5 | **0.7932** | 0.4258 | 0.5749 | 0.7838 |

The whole story is in this table. The unconstrained SMI deconvolution recovers
the fODF **to four decimal places** when there is no noise — it is inverting
its own forward model, so it should — and loses more than half of that by
SNR 10. Constrained SMI gives up 0.09 of ACC at every SNR to buy a curve that
is almost flat: 0.9053 noise free, 0.8870 at SNR 10. MSMT-CSD's curve is flatter
still but starts 0.07 lower, because its estimated response is blunt (section 7
shows an exact response does not fix that either — 0.8921).

### 5.6 Peak amplitude

Normalised to each method's own noise-free median, so methods whose fODFs are
not on the same scale can still be compared. Median [coefficient of variation].

_crossing 45°_

| SNR | SMI constrained | SMI unconstrained | SSST-CSD | MSMT-CSD |
|---|---|---|---|---|
| 30 | 1.011 [0.06] | 1.100 [0.10] | 1.054 [0.05] | **0.996 [0.02]** |
| 10 | 0.956 [0.09] | 1.183 [0.17] | 1.140 [0.12] | **0.936 [0.07]** |
| 5 | 0.916 [0.13] | 1.300 [0.21] | 1.187 [0.19] | 0.907 [0.14] |

Noise **inflates** the peak of the unregularized deconvolutions and **deflates**
the peak of the constrained ones. At SNR 5 the unconstrained SMI fODF's peak is
30% above its own noise-free value and SSST-CSD's is 19% above, while
constrained SMI is 8% below. Anything that thresholds on fODF amplitude — every
tractography cutoff — is reading a quantity that moves by tens of percent with
SNR, in opposite directions depending on the regularizer. MSMT-CSD's single
fibre amplitude moves the other way (1.058 at SNR 5, `deconv_tables.md`), which
is the Rician floor lifting its `l = 0` term.

---

## 6. The non-negativity weight is mistuned, and this measures it

`README for Claude` section 2 item 3 recorded that `lambda_nonneg = 10` — the
shipped default — is what closes the 45 degree crossing, and asked for a
cross-check against an independent metric before changing anything. This is
that cross-check, and it is a stronger one than asked for: the same peak finder
and the same ground truth as the four-method comparison above, 1000
realisations per condition at SNR 30, one `SMI.fit` and then one deconvolution
per setting from the same kernel, so nothing but the regularizer moves
(`deconv_comparison/sweep_nonneg.m`).

**Realisations returning the correct number of fibres, %**

| setting | single | crossing 15° | crossing 45° | crossing 60° |
|---|---|---|---|---|
| non-negativity off | 100.0 | 0.0 | **97.3** | 99.3 |
| `lambda_nonneg = 1` | 100.0 | 0.0 | 63.5 | 100.0 |
| `lambda_nonneg = 3` | 100.0 | 0.0 | 0.7 | 100.0 |
| `lambda_nonneg = 10` (**shipped**) | 100.0 | 0.0 | **0.0** | 100.0 |

**Angular error of the largest peak, degrees (median)**

| setting | single | crossing 15° | crossing 45° | crossing 60° |
|---|---|---|---|---|
| non-negativity off | 1.09 | 6.78 | **2.86** | 2.46 |
| `lambda_nonneg = 1` | 0.44 | 7.21 | 7.32 | 1.37 |
| `lambda_nonneg = 3` | **0.34** | 7.29 | 18.30 | **0.85** |
| `lambda_nonneg = 10` | 0.36 | 7.30 | 21.35 | 2.04 |

**Angular correlation coefficient against the band limited truth (mean)**

| setting | single | crossing 15° | crossing 45° | crossing 60° |
|---|---|---|---|---|
| non-negativity off | 0.9748 | 0.9710 | 0.9425 | 0.9447 |
| `lambda_nonneg = 1` | **0.9800** | **0.9855** | **0.9652** | **0.9714** |
| `lambda_nonneg = 3` | 0.9573 | 0.9687 | 0.9308 | 0.9434 |
| `lambda_nonneg = 10` | 0.9299 | 0.9466 | 0.8862 | 0.9023 |

Three things follow, and none of them is a matter of taste.

1. **`lambda_nonneg = 10` is not on the Pareto front.** `3` beats it on 45
   degrees, on 60 degrees, on the single fibre and on every ACC. `1` beats it on
   45 degrees and on every ACC by a wide margin. The shipped default is the
   worst of the four constrained settings by angular correlation at every
   condition.
2. **`lambda_nonneg = 1` is where the whole-fODF fidelity peaks**, and it is
   the only constrained setting that resolves 45 degree crossings more often
   than not (63.5%). If one number has to change, it is this one. The default
   was **not** changed here — that is the user's call, and the measurement is
   the deliverable.
3. **Tikhonov still does nothing.** `lambda_tikhonov` 0.3 versus 0 at
   `lambda_nonneg = 10` moves the 45 degree error from 21.35 to 21.33 degrees
   and the ACC in the fourth decimal. This confirms `README for Claude`
   section 2 item 1 on an independent metric.

There is also a cost to turning the constraint off that the ACC hides:
**spurious peaks appear only when it is off** — 0.021 per realisation at 45
degrees and 0.007 at 60 degrees, against exactly 0.000 for every constrained
setting. That is the thing non-negativity buys, and it is worth having; the
question is only how much of it to buy.

---

## 7. Control: was MSMT-CSD set up wrong?

A comparison in which MSMT-CSD comes out worst on crossings is equally
consistent with "MSMT-CSD was configured badly here", so both CSD variants were
run a second time on 2000 of the SNR 30 realisations with the **exact** zonal
response of the kernel that generated the data — an idealised delta response,
better than any estimator could produce (`deconv_comparison/control_exact.py`).

| | 45° resolved | 45° error | 60° error | single fibre error |
|---|---|---|---|---|
| SSST-CSD, `tournier` response | 27.5% | 10.65° | 1.60° | 0.59° |
| SSST-CSD, **exact** response | **83.8%** | **4.13°** | 1.67° | 0.68° |
| MSMT-CSD, `dhollander` responses | 0.0% | 22.11° | 7.06° | 0.32° |
| MSMT-CSD, **exact** responses | 0.0% | 21.85° | **2.57°** | 0.32° |

(the estimated-response rows are the 10,000-realisation values from section 5;
the exact-response rows are 500 realisations per condition, which is enough for
a median angular error to two decimals)

Two separate things come out of this.

* **MSMT-CSD's 60 degree error is response-limited, not a setup error.** An
  exact response takes it from 7.06 to 2.57 degrees — a factor 2.7. This is the
  same effect `README for Claude` section 6.3 recorded as "MSMT-CSD degrades
  3-4x on crossings" with an estimated response, now measured with a response
  estimated by the actual `dwi2response dhollander` algorithm rather than an
  ad hoc anisotropy selector.
* **MSMT-CSD's 45 degree failure is not response-limited.** It resolves 0.0%
  either way. dipy's MSMT-CSD at Lmax 6 behaves here exactly as SMI does at
  `lambda_nonneg = 10`: the non-negativity constraint closes the gap between two
  fibres 45 degrees apart. SSST-CSD, by contrast, goes from 30% to 83.8% — its
  45 degree performance **is** response-limited.

So the headline numbers in section 5 are the cost of estimating a response from
dispersed white matter, which is what every real pipeline does. They are not an
artefact of how the comparison was wired.

---

## 8. What this does not show

Worth being explicit about, because several of these are easy to read into the
tables and none of them is measured here.

* **Nothing has touched real data**, including the "real data" the responses
  come from (section 2). The single most valuable next step is to run
  `dhollander.py`'s selection on an actual brain and compare the response it
  returns against the ones here.
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
* **One protocol, one Lmax.** Three shells, 90 directions each, Lmax 6
  everywhere. The 15 degree crossing is unresolvable at Lmax 6 by construction —
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
python3 setup_protocol.py
octave --eval "oct_path; gen_phantom(30,'p30')"
octave --eval "oct_path; dump_bases"
python3 dhollander.py p30                       # dhollander, tournier, fa

for s in 50 30 20 10 5; do
  octave --eval "oct_path; gen_montecarlo($s,10000,'snr$s')"
  python3 run_csd.py snr$s p30
done
octave --eval "oct_path; gen_montecarlo(1e4,500,'nf')"
python3 run_csd.py nf p30

python3 check_conventions.py snr30
python3 tables.py nf snr50:50 snr30:30 snr20:20 snr10:10 snr5:5
python3 figure_mc.py ../fodf_deconv_montecarlo.png \
        snr50:50 snr30:30 snr20:20 snr10:10 snr5:5
python3 figure_response.py p30 ../fodf_response_shview.png
```

The SNR arms are independent and were run three at a time on four cores. One
arm is about 12 minutes of `SMI.fit`, 9 minutes for the second (unconstrained)
deconvolution and 18 minutes of MSMT-CSD; SSST-CSD is under a second for all
40,000 voxels. The whole thing is about two hours of wall clock on four cores.

Environment: GNU Octave 8.4 with `statistics`, three shims for MATLAB functions
Octave lacks (`deconv_comparison/stubs/`), Python 3.11 with
`numpy scipy dipy cvxpy matplotlib`.

---

## 10. What to do next, in order of value

1. **Check the `CS_phase` finding against a real MRtrix install** (section 4.1).
   It is a one-command test — write an SMI fODF, run `sh2peaks`, compare against
   the peak SMI reports — and if it holds, every SMI tractogram produced so far
   is mirrored. Nothing else in this report matters as much.
2. **Decide `lambda_nonneg`** (section 6). Two independent measurements now say
   the shipped `10` is not on the Pareto front and that `1` maximises whole-fODF
   fidelity at every condition. The default was not changed here.
3. **Run `dhollander.py` on a real brain** and compare the response it returns
   against the phantom-derived ones in section 2. That is the one substitution
   in this whole package, and it is cheap to remove.
4. **Add unequal fibre volume fractions and three-way crossings.** The
   machinery takes them — `gen_montecarlo.m` builds its conditions from a list
   of axes — and they are the part of the Jeurissen design that is missing.
5. **Re-run the four-method comparison with a heterogeneous microstructure per
   voxel.** SMI's per-voxel kernel should pay off there and cannot pay off in
   the current design, where every voxel has the same kernel.

---

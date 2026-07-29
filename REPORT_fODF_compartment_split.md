# Edema that spares the axons, response estimation, and the regularization sweeps

Measurement report. Follow-on to `REPORT_fODF_freewater.md`, answering three
questions: what happens when edema **redistributes** water instead of diluting
everything; what happens when CSD and MSMT-CSD have to **estimate** their
response functions from the most anisotropic voxels of a brain that contains
edema; and where SMI's angular resolution actually goes.

Simulation is `freewater_comparison/`. One of the three results contradicts a
conclusion in the previous report; that is called out in section 5.

---

## 1. The edema model

The previous report added free water by scaling everything down, so `f` fell
from 0.60 to 0.36 along with it. That is dilution. A more plausible account of
vasogenic edema is that the axons survive and the **extracellular** water is
what changes character: the intra-axonal fraction stays put and hindered
extra-axonal water becomes freely diffusing.

So `f = 0.60` is held **fixed** and the remaining 0.40 is walked across:

| | f | f_extra | fw |
|---|---|---|---|
| healthy | 0.60 | 0.40 | 0.00 |
| | 0.60 | 0.30 | 0.10 |
| **the case asked for** | **0.60** | **0.20** | **0.20** |
| | 0.60 | 0.10 | 0.30 |
| all extracellular water freed | 0.60 | 0.00 | 0.40 |

Protocol, noise, angular order, peak extraction and fairness controls are
unchanged from the previous report. Geometries are a single fibre and a 60
degree crossing, plus a CSF voxel, 40 realisations each.

**This distinction matters more than the free water fraction does.** Under
dilution, 40% free water cost CSD and MSMT-CSD 40% of their fODF amplitude.
Under redistribution, the same 40% free water costs CSD 13% and MSMT-CSD 23%:

| fw (f = 0.60 fixed) | SMI | CSD (b=3) | MSMT-CSD |
|---|---|---|---|
| 0.10 | 0.992 | 0.954 | 0.956 |
| **0.20** | **0.983** | **0.926** | **0.901** |
| 0.30 | 0.980 | 0.899 | 0.835 |
| 0.40 | 0.980 | 0.870 | 0.774 |

(median peak amplitude relative to the fw = 0 voxel of the same method,
single fibre, SNR 30; the 60 degree crossing behaves the same to within 0.02.)

The reason is that an apparent-fibre-density fODF tracks the **fibre** content,
not the total tissue fraction, and `f` is what did not move. SMI remains flat at
0.98-1.01 for the structural reason established previously: its fODF is
normalised to unit mass and carries no density at all.

**So the tracking problem in edema is much milder than the first report implied,
provided the axons survive.** How much milder depends entirely on how much of
the real perturbation is axonal loss versus water redistribution, which is a
question only real data can answer.

---

## 2. Response estimation from the most anisotropic voxels

CSD and MSMT-CSD do not get their response handed to them, as they did in the
previous report. They take it from the most anisotropic voxels of the volume.
To test whether that is safe in the presence of edema, two synthetic 4000-voxel
brains were built, identical except that one replaces 15% of its white matter
with the `f = 0.60, fw = 0.20` edema voxel above.

Two selectors, both taking the top 5%:

| selector | what it is |
|---|---|
| `FA` | tensor FA fitted on b <= 1. What dipy's `auto_response_ssst` and `mask_for_response_msmt` do. |
| `anisotropy` | `\|\|S_2\|\| / S_0` of the b = 3 shell signal itself. The high-b analogue. |

Composition of the selection:

| selector, brain | WM single | WM crossing | **EDEMA** | other |
|---|---|---|---|---|
| FA, healthy | 83.5% | 16.5% | – | 0.0% |
| FA, edema | 86.5% | 13.5% | **0.0%** | 0.0% |
| anisotropy, healthy | 83.0% | 16.5% | – | 0.5% |
| anisotropy, edema | 68.5% | 14.0% | **17.0%** | 0.5% |

(edema is 15% of that volume, so 17.0% is a mild enrichment, not exclusion.)

**On a high-b anisotropy criterion, edema-with-intact-axons is not
distinguishable from healthy white matter.** Its b = 3 shell anisotropy is
0.677 against healthy WM's 0.682; at b = 1 it is 0.396 against 0.395. Free water
is essentially gone by b = 3, so what remains is stick-plus-zeppelin at a
*higher* stick fraction (0.60/0.80 rather than 0.60/1.00), and the two effects
very nearly cancel. FA excludes it, but only through the low-b free water
signature, and by a margin (0.786 against 0.856) that a milder edema would close.

**And it does not matter.** The estimated response is the same either way:

```
WM response at b=3, normalised to l=0
  exact (delta fODF)     1.000  -0.828   0.441  -0.179
  FA,  healthy brain     1.000  -0.690   0.331  -0.098
  FA,  edema brain       1.000  -0.690   0.331  -0.096
  aniso, healthy brain   1.000  -0.707   0.338  -0.099
  aniso, edema brain     1.000  -0.707   0.341  -0.101
```

The contamination is benign precisely because the contaminant looks like the
thing being estimated. **What is not benign is the difference between any
estimated response and the idealised one**: every estimated response is far
blunter, because it absorbs the fibre dispersion of the voxels it was built
from. Handing CSD a delta-fODF response, as the previous report did, was not the
generous choice it was described as — it was a different operating point.

### What the bluntness costs

60 degree crossing, median peak angular error, SNR 30:

| fw | SMI | CSD (b=3), estimated | MSMT, exact | MSMT, estimated |
|---|---|---|---|---|
| 0.00 | **1.86** | 2.25 | 2.26 | 7.04 |
| 0.20 | **1.52** | 2.32 | 2.34 | 6.99 |
| 0.40 | **1.52** | 2.32 | 2.77 | 9.60 |

**CSD absorbs it; MSMT-CSD does not.** Single-shell CSD at b = 3 is unchanged
(2.2 vs 2.1 degrees). MSMT-CSD degrades 3-4x, and worsens as free water rises.
A per-shell tensor response (dipy's `response_from_mask_msmt` route) and a
per-shell per-l SH response (what MRtrix stores) were both tried and both show
it, 7.0 and 11.2 degrees respectively at fw = 0, so it is not dipy's tensor step.
The blunter response has smaller high-`l` coefficients, the deconvolution gain
at l = 6 roughly doubles, and MSMT's three-compartment QP handles that
amplification worse than CSD's simpler reweighting.

SMI has no response-estimation step to get wrong — it fits a kernel per voxel —
and is the most accurate of the four on every row.

MSMT with an estimated response does buy something: its CSF peak falls to 0.0026
and the edema/CSF contrast rises to 283x, against 24.7x with the exact response.
Better termination, worse orientation.

---

## 3. Compartment recovery, and why `f` is the interesting weight

Median, SNR 30. True `f` is 0.60 in every row but the last.

| f_extra / fw | SMI `f` | SMI `fw` | MSMT wm vf | MSMT csf vf |
|---|---|---|---|---|
| 0.40 / 0.00 | 0.515 | 0.023 | 1.060 | 0.000 |
| 0.30 / 0.10 | 0.596 | 0.082 | 1.016 | 0.000 |
| **0.20 / 0.20** | **0.638** | **0.151** | 0.966 | **0.019** |
| 0.10 / 0.30 | 0.633 | 0.270 | 0.892 | 0.092 |
| 0.00 / 0.40 | 0.617 | 0.384 | 0.818 | 0.166 |
| CSF | 0.114 | 1.000 | 0.004 | 0.995 |

Two things to read off this.

**MSMT-CSD's CSF volume fraction is not a usable free water estimate under this
model**: 0.019 where the truth is 0.20, and 0.166 where the truth is 0.40. With
`f_extra` still present the three-compartment model can explain the signal with
white matter alone. SMI's `fw` (0.151 at a true 0.20) is much closer, though it
is biased low throughout, partly because the fit absorbs free water into a
raised `Deperp` instead (0.48 -> 0.63 across the trajectory, noise-free).

**SMI's `f` is flat-to-rising across the whole edema trajectory.** That is the
property a modulation weight needs and that no `p2`-derived weight has.

### Weight comparison

Cutoff retains 95% of all WM including every free water level; the weight is
applied to SMI's fODF as a density modulation.

| weight | CSF surviving | edema (fw=.2) / healthy | fw=.4 / healthy |
|---|---|---|---|
| none | 5.0% | 0.98 | 0.98 |
| SMI `1-fw` | **0.0%** | 0.84 | 0.62 |
| **SMI `f`** | **0.0%** | **1.21** | **1.17** |
| SMI `p2product` | 0.0% | 0.94 | 0.92 |
| MSMT wm vol. frac. | **0.0%** | 0.95 | 0.80 |

**`f` suppresses CSF completely while leaving edema brighter than healthy white
matter.** Be careful about why: the true `f` is 0.60 everywhere, so the ideal
ratio is 1.00, and the 1.21 is an estimation bias (SMI reads 0.515 in the
healthy voxel and 0.638 in the edematous one), not a physical brightening. The
honest claim is that `f` is *flat* here, not that it favours edema.

This does **not** overturn `REPORT_fODF_modulation.md` section 7.5. That section
rejected `f` against an edema class with `f` collapsed to 0.30, where an
`f`-weight genuinely does delete the region of interest. The two results
bracket the question rather than contradict it:

> **A tissue-fraction weight fails on axonal loss and is safe under water
> redistribution.** Which regime real peritumoral tissue is in decides whether
> `f` is the right weight, and that is measurable in a known edema ROI.

---

## 4. The regularization sweeps

The previous report found that SMI never resolves a 45 degree crossing and
suspected `lambda_tikhonov`. Both regularizers were swept on identical data.

### Tikhonov does nothing to angular resolution

Noise-free, `lambda_nonneg = 10` throughout:

| lambda_tikhonov | 45° resolved | 60° error | CSF peak |
|---|---|---|---|
| 0.00 | 0/40 | 1.24° | 0.1015 |
| 0.05 | 0/40 | 1.24° | 0.1015 |
| 0.30 (shipped) | 0/40 | 1.24° | 0.1012 |
| 0.80 | 0/40 | 1.24° | 0.0998 |

Flat. The 45 degree crossing is lost at `lambda_tikhonov = 0`, with no noise.

### The non-negativity weight is the cause

Noise-free, `lambda_tikhonov = 0.3` throughout:

| lambda_nonneg | 45° resolved | 45° error | 60° error | CSF peak | WM/CSF |
|---|---|---|---|---|---|
| **off** | **40/40** | **1.69°** | 1.78° | 0.1625 | 8.77 |
| 1 | **40/40** | 6.52° | 1.78° | 0.1012 | **12.04** |
| 3 | 0/40 | – | 1.78° | 0.1012 | 10.35 |
| **10 (shipped)** | **0/40** | – | 1.24° | 0.1012 | 8.56 |
| 30 | 0/40 | – | 2.78° | 0.1012 | 7.75 |

With the constraint off, SMI resolves the 45 degree crossing in **40 of 40**
realisations at **1.69 degrees** — exactly the band-limited ground truth's own
error. At `lambda_nonneg = 1` it still resolves 40/40, and the WM/CSF contrast
is the best in the sweep. From 3 upward it never resolves.

**`lambda_nonneg = 10` is what closes the 45 degree crossing.** This does not
contradict `REPORT_fODF_regularization_sweep.md`, which measured 10 as the
optimum — that sweep optimised reconstruction error, not angular resolution, and
nobody had traded the two against each other. `lambda_nonneg = 1` looks like a
better operating point on this evidence, but the noise-free result alone is not
enough to recommend it: the non-negativity constraint is what prevents the CSF
blow-ups documented previously, and the SNR 30 and 15 arms of this sweep had not
finished at the time of writing. **Do not change the default on this table
alone.**

### Where the high-l power actually goes

Single fibre, no noise, rotational invariants against ground truth:

| | p2 | p4 | p6 |
|---|---|---|---|
| ground truth | 0.903 | 0.713 | 0.495 |
| SMI, lambda_tikhonov = 0 | 0.819 | 0.439 | **0.137** |
| SMI, lambda_tikhonov = 0.3 | 0.818 | 0.437 | **0.136** |

p6 comes back at 28% of truth and p4 at 62%, **with no noise and no Tikhonov
damping**, and the two lambdas are identical to three decimals. So
`REPORT_fODF_modulation.md` section 3's attribution of the `pl4` loss to the
Tikhonov term is **wrong**, at least at Lmax 6 and these lambdas. The loss is
monotone in `l` and is caused by the non-negativity constraint together with
error in the estimated kernel, whose `K_l` at high `l` is small and very
sensitive to the kernel parameters.

---

## 5. Corrections to the previous report

- The 45 degree resolution failure was attributed to `lambda_tikhonov`. It is
  `lambda_nonneg`. Tikhonov has no measurable effect on angular resolution at
  any value swept.
- CSD and MSMT-CSD were given delta-fODF responses and this was described as
  generous to them. It is better described as idealised: a realistically
  estimated response is much blunter, and for MSMT-CSD that costs 3-4x in
  crossing angular error.

## 6. Caveats

- Simulation only, on SMI's own forward model.
- The synthetic brains are scrambled voxel populations with no spatial
  structure, so they exercise the response estimators' *selection statistics*
  but not their ROI or spatial-contiguity heuristics.
- The `dwi2response tournier` peak-ratio criterion was implemented
  (`tournier.py`) but is not used for the headline numbers: at Lmax 6 a 50
  degree crossing does not resolve, presents a single peak, scores a
  second-to-first ratio of exactly 0, and is selected as the most
  single-fibre-like voxel in the volume — filling 88% of the selection with
  crossings in the healthy brain. That degeneracy is a real property of the
  criterion at this angular order and worth knowing, but it makes the resulting
  response useless as a comparison point.
- The SNR 30 and 15 arms of the non-negativity sweep were still running when
  this was written; only the noise-free arm is reported.

## 7. Reproducing

From `freewater_comparison/`, after `setup_protocol.py`:

```
octave-cli --no-gui run_c.m brain 30 healthy healthy30
octave-cli --no-gui run_c.m brain 30 edema   edema30
octave-cli --no-gui run_c.m cond  30 c_snr30       # and 15, 100000
python3 responses.py healthy30 && python3 responses.py edema30
python3 run_dipy2.py c_snr30                       # and c_snr15
python3 score2.py c_snr30 30
octave-cli --no-gui run_sweep.m clean              # lambda_tikhonov
octave-cli --no-gui run_nn.m    clean              # lambda_nonneg
python3 score_sweep.py && python3 score_nn.py
python3 figure2.py
```

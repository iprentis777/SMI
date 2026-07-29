# SMI vs CSD vs MSMT-CSD: recovering a WM fODF in 40% free water

Measurement report. Everything below is measured, not assumed; the simulation is
`freewater_comparison/`, and it is reproducible end to end from a bare container.

---

## 1. The question

How does the **regularized SMI deconvolution with no modulation** recover a white
matter fODF when the voxel is 40% free water, and how does that compare with
single-shell CSD and MSMT-CSD?

40% free water with intact fibres is the peritumoral edema case in its mildest
honest form: the tissue microstructure is held **fixed** and water is added, so
this isolates the free water effect from the axonal loss that usually accompanies
it. It is deliberately a weaker perturbation than the edema class in
`REPORT_fODF_modulation.md`, which also collapses `f` from 0.60 to 0.30.

---

## 2. Design, and why it is fair

| | |
|---|---|
| Protocol | HCP-like, b = 0/1/2/3 ms/um^2, 18/90/90/90 directions, 288 volumes |
| Noise | Rician, SNR 30 and 15 relative to b=0; 40 realisations per condition |
| Angular order | **Lmax 6 for all three methods** |
| Conditions | single fibre, 60 deg crossing, 45 deg crossing, each at fw = 0 and fw = 0.40; plus pure CSF |
| WM tissue | `Da 2.0, Depar 2.0, Deperp 0.50`, intra-axonal fraction 0.60 within tissue, Watson kappa 16 |

Three deliberate choices remove the usual confounds:

1. **One signal, three fitters.** The DWI is synthesised once by SMI's own
   forward model (`fODF_modulation_helpers.m`) and written to disk. Octave/SMI
   and Python/dipy read the *same bytes*. No re-implementation can drift.
2. **One peak finder, one sphere.** Every method's fODF is evaluated on the same
   11554-vertex sphere (1.95 deg spacing) and peaks are extracted by the same
   `dipy.direction.peaks.peak_directions` call. No result can come from a
   different peak finder.
3. **CSD and MSMT-CSD are given the exact ground-truth response functions**,
   derived analytically from the same kernel that generated the data
   (`kernel.py`, verified against `SMI.RotInv_Kell_wFW_b_beta_TE_numerical` to
   `2.6e-15`). They carry **zero response-estimation error**. SMI has to estimate
   its kernel from the noisy data. The comparison is therefore biased *against*
   SMI, on purpose.

The WM response given to CSD and MSMT is the **fw = 0** tissue response — what
`dwi2response tournier` would return from healthy single-fibre white matter.
Neither is told that the voxel it is handed contains 40% free water.

**Validation.** With no noise, MSMT-CSD recovers volume fractions of
WM 0.628 / CSF 0.368 / GM 0.000 against a true 0.60 / 0.40 / 0, and every method
lands within 0.75 deg of the ground truth for a single fibre. The construction is
sound.

---

## 3. The headline result

**SMI's fODF amplitude is essentially blind to free water. CSD's and MSMT-CSD's
scale with the tissue fraction.**

Peak amplitude ratio, `fw = 0.40` over `fw = 0`, SNR 30:

| geometry | **SMI** | CSD (b=3) | CSD (b=1) | MSMT-CSD |
|---|---|---|---|---|
| single fibre | **0.96** | 0.60 | 0.62 | 0.62 |
| 60 deg crossing | **1.01** | 0.62 | 0.65 | 0.61 |
| 45 deg crossing | **0.94** | 0.62 | 0.62 | 0.59 |

CSD and MSMT-CSD all land on `1 - fw = 0.60`. That is not a defect — it is what
an apparent-fibre-density fODF is *supposed* to do. SMI lands on 1.0 because
`out.plm` is stored in the normalized convention `p_00 = 1`, so the fODF
integrates to 1 in every voxel and carries no density information to lose.

At SNR 15 the same holds (SMI 0.94 / 1.00 / 0.86; the others 0.59-0.69).

### What this costs, and what it buys

Median peak fODF amplitude, SNR 30:

| condition | **SMI** | CSD (b=3) | MSMT-CSD |
|---|---|---|---|
| healthy WM | 0.855 | 1.255 | 0.787 |
| WM + 40% free water | **0.824** | 0.750 | 0.486 |
| CSF | **0.321** | 0.076 | 0.028 |

| SNR 30 | **SMI** | CSD (b=3) | CSD (b=1) | MSMT-CSD |
|---|---|---|---|---|
| (WM + 40% FW) / CSF peak ratio | 2.57 | 9.82 | 9.52 | **17.14** |
| CSF above MRtrix's default 0.05 | 100% | 100% | 90% | **0%** |
| CSF above a cutoff keeping 95% of WM | 5.0% | 0% | 0% | **0%** |

| SNR 15 | **SMI** | CSD (b=3) | CSD (b=1) | MSMT-CSD |
|---|---|---|---|---|
| (WM + 40% FW) / CSF peak ratio | 1.99 | 4.96 | 6.40 | **8.89** |
| CSF above MRtrix's default 0.05 | 100% | 100% | 100% | 87.5% |
| CSF above a cutoff keeping 95% of WM | 20.0% | 0% | 0% | **0%** |

So the two properties are the same property seen from opposite ends:

- **SMI is free-water-invariant.** A voxel with 40% free water tracks exactly as
  well as a healthy one. But CSF is at 0.321 — four times the `1/(4*pi) = 0.0796`
  isotropic floor, because the deconvolution amplifies noise where `K_2 -> 0`,
  even regularized — and 100% of CSF survives MRtrix's default cutoff.
- **MSMT-CSD is free-water-proportional.** CSF is annihilated (0.028, 0% above
  cutoff), but the WM fODF in a 40%-free-water voxel is dimmed to 0.62 of
  healthy, so a fixed cutoff tuned on healthy WM will preferentially terminate
  tracks in edema.

This is the trade-off the whole modulation effort has been circling, now measured.
The unmodulated SMI fODF is the *only* one of the three that does not dim in
edema, and it is also the only one that cannot terminate in CSF.

---

## 4. Orientation accuracy

Median peak angular error in degrees (IQR), SNR 15. `nres` counts realisations
that failed to resolve all fibres.

| condition | **SMI** | CSD (b=3) | CSD (b=1) | MSMT-CSD |
|---|---|---|---|---|
| single, no FW | **0.75** (0.55) | 1.30 (0.99) | 1.65 (1.12) | **0.75** (0.55) |
| single, 40% FW | **1.30** (0.99) | 2.60 (1.89) | 2.61 (2.49) | **1.30** (1.22) |
| 60 deg cross, no FW | **2.71** (1.82) | 3.28 (2.19) | 5.51 (3.42) | 2.68 (1.86) |
| 60 deg cross, 40% FW | **3.50** (2.38) | 5.64 (3.63) | 9.70 (4.07) | 3.51 (2.63) |
| 45 deg cross, no FW | *unresolved* 40/40 | 4.86 (3.00), 16 nres | 7.32, 33 nres | *unresolved* 40/40 |
| 45 deg cross, 40% FW | *unresolved* 40/40 | 8.95 (8.13), 12 nres | 8.27, 31 nres | *unresolved* 40/40 |

Two findings:

1. **On orientation, SMI and MSMT-CSD are indistinguishable and both beat CSD**,
   at both free water levels. Adding 40% free water roughly doubles the angular
   error for every method — that part is unavoidable, it is a real loss of SNR in
   the anisotropic signal.
2. **SMI cannot resolve a 45 degree crossing, in any realisation, at any SNR —
   not even with no noise at all.** The band-limited ground truth resolves it
   (peaks 1.69 deg from truth), and single-shell CSD at b=3 resolves it in 24/40
   realisations at SNR 30. So this is not an information limit, it is the
   regularized deconvolution's angular resolution. The likely cause is
   `lambda_tikhonov = 0.3` damping the l=4 and l=6 bands, which is the same
   mechanism `REPORT_fODF_modulation.md` section 3 identified when it found `pl4`
   was measuring the regularizer rather than the tissue. MSMT-CSD also fails at
   45 deg, for a different reason: it averages over all four shells, including
   b=1 where angular contrast is weakest.

---

## 5. Free water estimation

Both frameworks estimate the free water fraction well; they differ only in what
they then do with it.

| condition (true fw) | SMI `fw` | MSMT CSF vf | MSMT GM vf | MSMT WM vf |
|---|---|---|---|---|
| single, fw = 0 | 0.018 | 0.000 | 0.000 | 1.034 |
| single, fw = 0.40 | **0.365** | **0.362** | 0.000 | 0.634 |
| 60 deg cross, fw = 0.40 | 0.347 | 0.380 | 0.000 | 0.619 |
| CSF (fw = 0.95) | 1.000 | 0.909 | 0.000 | 0.087 |

(SNR 30. GM is a nuisance compartment — the ground truth contains none, and MSMT
correctly assigns it zero, so there is no leakage.)

---

## 6. An unexpected result on tissue-fraction weights

Because both methods estimate free water accurately, it is worth asking directly
whether a tissue-fraction weight works *in this specific condition*. Cutoff chosen
to retain 95% of all WM including the 40%-free-water voxels, then CSF survival:

| weight | CSF surviving, SNR 30 | SNR 15 | median edema / healthy amplitude |
|---|---|---|---|
| none (unmodulated) | 5.0% | 20.0% | 0.87-0.93 |
| SMI `1-fw` | **0.0%** | **0.0%** | 0.60 |
| SMI `f` | **0.0%** | **0.0%** | 0.71-0.74 |
| SMI `p2product` (shipped default) | 0.0% | 7.5% | 0.70-0.88 |
| MSMT WM volume fraction | **0.0%** | **0.0%** | 0.55-0.57 |

**This does not overturn `REPORT_fODF_modulation.md` section 7.5, and should not
be read as doing so.** That section rejected `f` and `1-fw` against an edema
class with `f` collapsed to 0.30 *and* `fw = 0.50`, where a tissue-fraction weight
genuinely does delete the region of interest. Here the tissue is intact and only
diluted, and a tissue-fraction weight merely dims it by ~40% rather than removing
it. The two results measure different severities and are both correct.

What it does establish is the **boundary** of the objection: tissue-fraction
weights fail on axonal loss, not on free water as such. Distinguishing those two
is what a real-data experiment in a known edema ROI would settle.

---

## 7. What this implies

1. **Nothing here argues for switching away from SMI on orientation grounds.**
   SMI matches MSMT-CSD and beats CSD on peak angular error, while estimating its
   own kernel rather than being handed one.
2. **The CSF problem is entirely an amplitude-convention problem, and it is
   SMI-specific.** MSMT-CSD gets CSF suppression for free because its fODF
   carries density. SMI's does not, by construction. This is the strongest
   argument yet that the fix belongs at export — either a weight or the
   `A(u) - 1/(4*pi)` floor subtraction — rather than in the deconvolution.
3. **The 45 degree resolution failure is new and is not about free water.** It is
   present with zero noise and is worth a Tikhonov sweep at Lmax 6 on its own.
   `lambda_tikhonov` was tuned for stability, not for angular resolution, and no
   measurement has yet traded the two off against each other.
4. **`p2product` is again the weakest of the candidate weights at SNR 15** (7.5%
   CSF surviving against 0.0% for the tissue fractions), consistent with the
   handoff's judgement that it should not be treated as endorsed.

### Caveats

- Simulation only. Every number here assumes the standard model kernel is the
  truth, which is exactly the assumption SMI makes, so SMI is being scored partly
  on its home ground even though its kernel is estimated.
- The "edema" here is pure dilution with intact fibres. Real peritumoral tissue
  has axonal loss too, and section 6 above is sensitive to which of the two
  dominates.
- Single-shell CSD is run at b = 3 ms/um^2 (primary) and b = 1 (secondary).
  b = 1 is worse on every metric, as expected: free water is barely attenuated
  there.
- MSMT-CSD is run as standard 3-tissue WM/GM/CSF. dipy's `multi_tissue_basis`
  rejects a 2-compartment configuration, so an exactly-matched WM+FW model was
  not available.

---

## 8. Reproducing it

Needs Octave (with `statistics` and `image`) and Python with
`numpy scipy dipy cvxpy matplotlib`. From `freewater_comparison/`:

```
python3 setup_protocol.py                     # protocol + shared evaluation sphere
octave-cli --no-gui run_smi.m 100000 clean    # noise-free reference
octave-cli --no-gui run_smi.m 30 snr30
octave-cli --no-gui run_smi.m 15 snr15
python3 run_dipy.py clean && python3 run_dipy.py snr30 && python3 run_dipy.py snr15
python3 validate.py                           # noise-free sanity table
python3 score.py snr30 30 && python3 score.py snr15 15
python3 weights.py
python3 figure.py                             # fodf_freewater_comparison.png
```

`SMI.fit` takes about 60 s per SNR for the 280-voxel volume; the MSMT-CSD QP
dominates the Python side.

| file | contents |
|---|---|
| `setup_protocol.py` | gradient scheme and the shared 11554-vertex evaluation sphere |
| `gen_and_fit_smi.m` | ground truth, signal synthesis, the regularized SMI fit |
| `kernel.py` | the SM kernel `K_l(b)` in Python, for building responses only |
| `check_kernel.m` | verifies `kernel.py` against SMI's Octave original |
| `run_dipy.py` | CSD and MSMT-CSD with exact ground-truth responses |
| `common_score.py` | the single shared peak extractor |
| `validate.py` / `score.py` / `weights.py` | noise-free check, main tables, weight comparison |
| `figure.py` | `fodf_freewater_comparison.png` |
| `stubs/` | Octave shims for `round(x,n)`, `discretize`, `datetime` |

Octave notes, in addition to those in `README for Claude.md` section 4: the SM
kernel hard-codes Legendre polynomials only to `l = 8`
(`SMI.m:2099-2107`), which caps the order at which ground truth can be
synthesised. `SMI.fit` defaults to `compartments = {'IAS','EAS'}`, so
**`'FW'` must be requested explicitly** or free water is not fitted at all.

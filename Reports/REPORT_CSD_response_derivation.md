# The response function given to the CSD arms: derivation

How the response handed to `dwi2fod csd` and `dwi2fod msmt_csd` in
`deconv_comparison/notebooks/smi_manuscript_60deg.m` is derived, why it takes
the form it does, and what it assumes. Written at the level of detail a
manuscript methods section needs, because the response is the one object both
methods share and the comparison rests on it being right.

Everything here is verified numerically inside the simulation on every run; the
checks are named at the end.

---

## 1. The problem the derivation solves

SMI and CSD do not describe tissue with the same object.

- **SMI has a kernel.** The Standard Model compartment vector
  `[f, Da, Depar, Deperp, fw]`, whose response to diffusion weighting enters the
  forward model through rotational invariants `K_l(b)`.
- **CSD has a response function.** A table of zonal spherical harmonic
  coefficients, one row per shell and one column per even `l`, which MRtrix
  reads from a `.txt` file.

To compare the two on one dataset, one has to be expressed in the other's terms.
The conversion is exact, it is one line, and section 3 derives it.

---

## 2. Conventions

Both conventions are fixed and neither is negotiable; getting either wrong
produces a plausible-looking fODF that is wrong by a rotation or a scale.

### 2.1 SMI's forward model

For a voxel whose fibre orientation distribution has normalised coefficients
`p_lm` (SMI's convention: `p_00 = 1`, so the fODF integrates to 1 over the
sphere), the signal in direction `u` at b-value `b`, normalised by `S0`, is

```
S(u)/S0 = sum_{l,m} K_l(b) · p_lm · Y_lm(u) · sqrt((2l+1)·4π)          (1)
```

This is `SMI.m:818-820`, and it is what `SMI.get_plm_from_S_and_kernel`
inverts. `Y_lm` is the real, even-order spherical harmonic basis returned by
`SMI.get_even_SH`.

### 2.2 The kernel's rotational invariants

For linear tensor encoding at a single echo time, the Standard Model with a free
water compartment gives

```
K_l(b) = ∫₀¹ [ f·e^(−b·Da·x²)
             + (1−f−fw)·e^(−b·(Deperp + (Depar−Deperp)·x²))
             + fw·e^(−b·D_FW) ] · P_l(x) dx                            (2)
```

with `x = cos θ` the cosine of the angle to the fibre axis and `P_l` the
Legendre polynomial. This is `SMI.RotInv_Kell_wFW_b_beta_TE_numerical`
(`SMI.m:2415`), evaluated by 200-point Gauss–Legendre quadrature. There is no
`(2l+1)` prefactor. The free-water term is carried at every `l` and removes
itself for even `l ≥ 2`, because `∫₀¹ P_l(x) dx = 0` there — free water is
isotropic and contributes only to `l = 0`, as it must.

### 2.3 MRtrix's response convolution

MRtrix treats the signal as the fODF convolved with an axially symmetric
response. In spherical harmonics that convolution is a per-band product:

```
S_lm = f_lm · r_l · sqrt(4π/(2l+1))                                     (3)
```

where `f_lm` are the fODF's SH coefficients and `r_l` the response's zonal
coefficients — exactly the numbers stored in a response `.txt`, one row per
shell, one column per even `l`, with amplitude reconstructed as `Σ_l r_l Y_l0`.

**Equation (3) is the single most important structural fact in this report,**
and section 6 returns to it: convolution acts on each band `l` by a scalar. It
cannot mix `m` within a band.

---

## 3. The delta response

Take a single fibre with no dispersion — a delta fODF on the sphere, pointing
along `z`. In SMI's normalised convention every zonal coefficient of a delta is
`p_l0 = 1`. Substituting into (1) for that voxel,

```
S(θ)/S0 = sum_l K_l(b) · Y_l0(θ) · sqrt((2l+1)·4π)
```

and matching term by term against the response form `Σ_l r_l Y_l0(θ)` gives

```
r_l(b) = K_l(b) · sqrt((2l+1)·4π)                                       (4)
```

This is implemented once, in `helpers/SMI_response_helpers.m` (`H.zh`), and it
is the bridge between the two toolboxes. A response file written from (4) can be
handed to `shview`, `dwi2fod` or `amp2response` with no further conversion.

**Verification.** The simulation reconstructs the single-fibre signal twice —
once through (1) using SMI's own forward code, once through (4) as a response
profile — and requires agreement. Measured `max|err| = 1.03e-15`
(`CHECK zonal response == SMI forward model`, Step 3).

---

## 4. The dispersion correction

Equation (4) describes a fibre with no dispersion. **No fibre in this simulation
is one.** Each population is a Watson distribution with concentration
`κ = 16` (about 14° of dispersion), chosen because a response estimated from
real white matter has already absorbed fibre dispersion, and a delta ground
truth would create a response/truth mismatch that does not exist in practice.

The correct response is therefore the response of *one dispersed population*,
not of one delta. Because the ground truth is built as

```
fODF = Watson(a₁, κ) + Watson(a₂, κ)  =  Watson(·, κ) ⊛ (δ_{a₁} + δ_{a₂})
```

and convolution is associative, the signal can be regrouped as

```
signal = kernel ⊛ Watson ⊛ (δ_{a₁} + δ_{a₂})
       = [ the dispersed response ] ⊛ [ a two-delta fODF ]
```

By (3), convolving two axially symmetric functions multiplies their normalised
zonal coefficients. Writing `p_l^W(κ)` for the Watson's own normalised zonal
coefficients (`p_0^W = 1`),

```
r_l(b) = K_l(b) · p_l^W(κ) · sqrt((2l+1)·4π)                            (5)
```

This is `RESPONSE_MODE = 'dispersed'`, the shipped default. `'delta'` recovers
(4).

### 4.1 The Watson zonal coefficients

For an axially symmetric distribution `f(x)`, `x = cos θ`, the normalised zonal
coefficients in the `p_0 = 1` convention are

```
p_l^W = ∫₀¹ e^(κx²) · P_l(x) dx  /  ∫₀¹ e^(κx²) dx                      (6)
```

The `[−1, 0]` half adds nothing: `f` is even in `x` and so is every even `P_l`.
At `κ = 16`,

| `l` | 0 | 2 | 4 | 6 | 8 |
|---|---|---|---|---|---|
| `p_l^W` | 1 | 0.9027 | 0.7126 | 0.4945 | 0.3037 |

These are the factors that blunt the delta response. A value of 1 at every `l`
would be a delta.

**Verification.** The simulation computes (6) two independent ways — by
projecting a sampled Watson onto spherical harmonics on the quadrature grid, the
same code path the ground truth uses, and by direct 1-D Gauss–Legendre
quadrature with no spherical harmonics anywhere. They agree to
`max|err| = 3.4e-05`, which is the band-limiting error of the SH route and not a
disagreement (`CHECK Watson p_l, harmonics vs 1-D quadrature`, Step 2).

### 4.2 Where the dispersed response sits

Normalised to `l = 0` at `b = 3 ms/µm²`:

| response | `r_2` | `r_4` | `r_6` |
|---|---|---|---|
| delta, equation (4) | −0.828 | 0.441 | −0.179 |
| **dispersed, equation (5)** | **−0.748** | **0.315** | **−0.089** |
| `dwi2response dhollander` | −0.703 | 0.342 | −0.107 |
| `dwi2response tournier` | −0.714 | 0.349 | −0.116 |
| `dwi2response fa` | −0.718 | 0.352 | −0.121 |

At `l = 2` the dispersed response lands between the delta and the three
estimators, accounting for about **69%** of the gap between them. That is a
quantitative statement of something the literature asserts qualitatively: most
of what makes an estimated response blunter than an idealised kernel is fibre
dispersion in the voxels the response was estimated from. At `l = 4` and `l = 6`
it slightly overshoots, so `κ = 16` disperses a little more than the phantom
those three were measured on.

The estimator rows are indicative only — they were measured on a superseded
synthetic protocol, not the acquisition used here.

---

## 5. Which tissue's kernel, and why it is the healthy one

Equation (5) needs a kernel. The simulation contains two tissues:

| preset | `[f, Da, Depar, Deperp, fw]` |
|---|---|
| healthy | `[0.60, 2.0, 2.0, 0.50, 0.02]` |
| edema | `[0.10, 2.4, 2.7, 1.15, 0.35]` |

**Both are deconvolved with the response built from the healthy kernel**
(`CSD_RESPONSE_KERNEL = 'healthy'`).

This is a statement about workflows, not about kernels. A CSD response is
estimated once per subject or per study by selecting single-fibre white matter
voxels and averaging them — `dwi2response tournier` and `dhollander` both work
by finding the most anisotropic voxels available. Nobody estimates a response
*for* edema: the edema is the tissue under examination rather than a reference
population, and it is close to the least anisotropic tissue in the brain, so the
selection heuristics would exclude it by construction. On real data, edematous
voxels are deconvolved with a healthy-white-matter average, and the mismatch
between that response and the tissue actually present is an error the method
carries.

Giving the edema arm an edema-derived response would hand CSD information it
cannot have in practice, and would make the edema comparison flattering to it.
That configuration remains available as `CSD_RESPONSE_KERNEL = 'matched'`, and
is worth running as a control, because it separates *"CSD is hurt by the
response mismatch"* from *"CSD is hurt by the low anisotropic signal in edema"*.

**This is the one place the two arms are not on equal footing, and the
inequality is a real property of the methods rather than an artefact of the
simulation.** SMI estimates its kernel per voxel and therefore adapts to edema —
imperfectly, and Step 6 prints how imperfectly. CSD is handed a fixed response
and cannot adapt. That asymmetry *is* the comparison.

### 5.1 A consequence for MSMT-CSD specifically

MSMT-CSD separates tissues by how the `l = 0` response decays across shells. The
edema kernel's `l = 0` decay is close to that of free water:

| `l = 0`, normalised | b = 0 | b = 1 | b = 2 | b = 3 |
|---|---|---|---|---|
| healthy WM | 1 | 0.513 | 0.335 | 0.253 |
| edema WM | 1 | 0.188 | 0.069 | 0.040 |
| CSF | 1 | 0.050 | 0.003 | 0.000 |

The angle between the WM and CSF response vectors falls from **31.3°** (healthy)
to **8.9°** (edema). At 8.9° those two columns of MSMT's design matrix are
nearly parallel, the tissue split becomes ill-conditioned, and the recovered WM
fODF degrades: on noise-free edema signal, MSMT recovers a 60° crossing at
60.94° with the healthy response and at **21.49°** with the matched one.

So the realistic choice is also the better-conditioned one. This is a
non-obvious result and it runs opposite to the intuition that a
better-matched response must give a better answer.

---

## 6. What the response does and does not control

Equation (3) makes deconvolution a **per-band scalar division**:

```
f_lm = S_lm / ( r_l · sqrt(4π/(2l+1)) )
```

Using a wrong response `r'_l` in place of `r_l` therefore multiplies each band
by the scalar `r_l / r'_l`. **It cannot mix `m` within a band, so it cannot
rotate anything.** Three consequences, all measured:

1. **For a single fibre, peak orientation is invariant to the response.** The
   fODF is axially symmetric about the fibre axis, and any per-band rescaling
   leaves it axially symmetric about the *same* axis. Measured on noise-free
   edema signal, the recovered peak error is **1.721°** whether the response is
   matched, healthy, has its `l ≥ 2` bands halved, or has them tripled —
   identical to three decimals.

2. **For a crossing, orientation is near-invariant and the lobes are not.**
   Sweeping the `l ≥ 2` scaling over a 16-fold range on a 60° crossing with
   unequal fibre weights (0.70 / 0.30) moves the dominant peak by at most 1.1°,
   while the secondary peak moves by 2.6°, the apparent separation by 3°, and at
   the lowest scaling the second fibre is **not detected at all**.

3. **Amplitude is not invariant.** The healthy and edema responses differ by
   **6.3×** in absolute `l = 0` amplitude but only 4–27% in normalised shape.
   The fODF absorbs the amplitude difference, which is why orientation survives
   — and why amplitude is where the mismatch actually lands.

**The methodological consequence.** A comparison scored only on angular error
will find CSD almost immune to response mis-specification, and will therefore
miss the effects that matter for tractography: amplitude collapse, loss of the
weaker fibre, and unstable tissue fractions. The simulation scores fibre count,
spurious peaks and peak amplitude alongside angular error for this reason.

---

## 7. The isotropic responses for MSMT-CSD

`dwi2fod msmt_csd` requires at least as many tissue responses as there are
shells to separate. The simulation contains one tissue, so the second and third
are idealised isotropic responses:

```
r_0(b) = e^(−b·D) · sqrt(4π),     D = D_GM = 0.8 µm²/ms   (grey matter)
r_0(b) = e^(−b·D) · sqrt(4π),     D = D_FW = 3.0 µm²/ms   (free water)
```

One column each, hence `l = 0` only, which is also what fixes their angular
order at 0 without any `-lmax` argument. The `sqrt(4π)` converts an isotropic
amplitude `A` to its zonal coefficient, since `Y_00 = 1/sqrt(4π)`.

**No simulated voxel contains grey matter or free-standing CSF.** Whatever MSMT
assigns to those compartments is therefore a measurement of its own leakage, and
it is reported rather than discarded. `D_FW` is the same value the forward model
uses for its free-water compartment, so the CSF response is consistent with the
simulation's own physics by construction. `D_GM = 0.8 µm²/ms` is a conventional
literature value and is **the one constant here not derived from anything else
in the package** — no simulated voxel is grey matter, so it cannot bias the
white matter fODF much, but it is an input rather than a result.

---

## 8. Assumptions

Listed because a methods section has to, and because each hands a method
something it would not have on real data.

1. **The response is exact, not estimated.** Equation (5) is evaluated
   analytically from a known kernel. A real workflow runs `dwi2response` and
   inherits both estimation noise and the spread of the voxels it averaged over.
   This makes the CSD arms *better* than they would be in practice, in a way
   that is not quantified here.
2. **The dispersion in the response is the true `κ`.** A real estimate absorbs
   whatever dispersion the selected voxels happen to have.
3. **The response is evaluated at MRtrix's own per-shell b values**, read back
   with `mrinfo -shell_bvalues`, not at nominal values — this protocol's shells
   jitter over 18 distinct b values. But MRtrix deconvolves with one response
   per shell while the signal was generated at each volume's exact b; that
   residual is measured and printed in Step 6b.
4. **`D_GM` is assumed**, per section 7.
5. **The generating model and SMI's fitted model are the same Standard Model.**
   There is no model mismatch anywhere in this simulation. This is the largest
   assumption in the package and it favours SMI.

---

## 9. Where each claim is checked

| claim | check | measured |
|---|---|---|
| equation (4) is consistent with SMI's forward model | `CHECK zonal response == SMI forward model`, Step 3 | 1.03e-15 |
| equation (6) is right | `CHECK Watson p_l, harmonics vs 1-D quadrature`, Step 2 | 3.4e-05, band limiting |
| the response file holds the numbers computed | `CHECK Lmax n response round trip`, Step 6b | < 1e-7 |
| MRtrix and the simulation agree on the shells | `CHECK MRtrix and SMI agree on the shells`, Step 6b | exact |
| the response order is what was intended | column count of the response file; verified against `dwi2fod` output coefficients | 15 / 28 / 45 at Lmax 4 / 6 / 8 |
| the mismatch when the response is not the tissue's | printed in Step 6b as `MISMATCH` | per run |

`deconv_comparison/test_csd_arms.m` re-runs the end-to-end consequence in about
two seconds.

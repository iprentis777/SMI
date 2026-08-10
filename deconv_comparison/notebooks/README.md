# Notebooks

Narrated, step-by-step versions of the simulations, for reading rather than
running in bulk. They are the answer to "how do I check this myself?".

| notebook | what it walks through |
|---|---|
| `smi_simulation_walkthrough.m` | the SMI arm end to end on a real HCP protocol: ground truth, kernel, forward convolution, Rician noise, `SMI.fit` at Lmax 4/6/8, peaks, MRtrix export. Conditions: single fibre, 30, 45, 60 degrees. **One SNR** |
| `smi_manuscript_60deg.m` | **the manuscript figure source.** 60 degree crossing only, two kernels back to back, swept across SNR, and **all three arms — SMI, SSST-CSD and MSMT-CSD — run inside this one file** |
| `smi_manuscript_60deg.ipynb` | an earlier three-arm version as a Jupyter notebook on the Octave kernel. Healthy kernel only, exact delta response. **Superseded by the `.m` for the manuscript**; kept because its figures are verified to render inline |

## `smi_manuscript_60deg.m` runs all three arms itself

The `.m` file calls the MRtrix3 binaries directly, in Step 6b, on **the same
`S_noisy` array `SMI.fit` was just given**.

| arm | what it is |
|---|---|
| SMI | `SMI.fit`, `fODF_regularization.flag_nonneg = 1`, `CS_phase = 0` |
| SSST-CSD | `dwi2fod csd` on the top shell — the MRtrix3 binary |
| MSMT-CSD | `dwi2fod msmt_csd` on all four shells, 3 tissues — the binary |

That single fact is what the comparison rests on: **one simulation, one noise
draw, one peak finder, one file.** There is no second forward model to keep in
step, no second noise stream, and no question about whether the arms saw the
same data — they are the same array. `RUN_MRTRIX = 0` scores the SMI arm alone.

Nothing here reimplements MRtrix. `dwi2fod`, `dwiextract` and `mrinfo` are the
real binaries; the only MRtrix behaviour implemented locally is reading and
writing its image format (`mrtrix_io.m`), and `mrinfo` checks that on every run.

### The response the CSD arms are given is dispersion matched

SMI has no response function, it has a kernel; the conversion is one line. But
the *delta* response `r_l(b) = K_l(b) sqrt((2l+1) 4 pi)` describes a fibre with
no dispersion, and no fibre in this simulation is one — each population is a
Watson at `kappa = 16`. So `RESPONSE_MODE = 'dispersed'` (the default) builds

```
r_l(b) = K_l(b) * p_l^Watson(kappa) * sqrt((2l+1) * 4*pi)
```

the response of *one such population*, which is what a perfect `dwi2response`
run on this data would recover. `'delta'` switches back. Measured at b = 3,
normalised to `l = 0`:

| response | `r_2` | `r_4` | `r_6` |
|---|---|---|---|
| delta (exact kernel) | -0.828 | 0.441 | -0.179 |
| **dispersion matched, `kappa = 16`** | **-0.748** | **0.315** | **-0.089** |
| `dwi2response dhollander` | -0.703 | 0.342 | -0.107 |
| `dwi2response tournier` | -0.714 | 0.349 | -0.116 |
| `dwi2response fa` | -0.718 | 0.352 | -0.121 |

At `l = 2` it lands **between** the delta and the estimators, accounting for
about **69%** of the gap — most of what makes an estimated response blunter than
the kernel really is fibre dispersion, which is what section 6.3 of the Monte
Carlo report supposed. At `l = 4` and `l = 6` it overshoots slightly, so
`kappa = 16` disperses a little more than the phantom those estimators ran on.
Read that row as indicative: they were estimated on the **superseded synthetic
protocol**, not this acquisition.

Note that `RESPONSE_MODE` decides the *shape* of the response (delta vs
dispersed). **Which tissue's kernel it is built from is a separate setting**,
`CSD_RESPONSE_KERNEL`, and its default is not this iteration's kernel — see the
next section.

### Which tissue's kernel the CSD response comes from

`CSD_RESPONSE_KERNEL` decides this, and it is the assumption the manuscript is
really about. The default is **`'healthy'`: both tissues are deconvolved with
the healthy WM response.** `'matched'` gives each tissue a response built from
its own kernel, and is a control rather than a realistic configuration.

**Why `'healthy'` is the realistic setting.** A CSD response is estimated once
per subject or per study by selecting single-fibre white matter voxels and
averaging them. Nobody estimates a response *for edema*: it is the tissue being
imaged, not a reference population, and the selection heuristics
(`dwi2response tournier`, `dhollander`) look for the most anisotropic voxels,
which is precisely what edema is not. So on real data edema voxels are
deconvolved with a healthy-WM average, and the mismatch is an error the method
has to live with.

**This is the one place the two arms are not on equal footing, and the
inequality is real rather than an artefact of the setup.** SMI re-estimates its
kernel *per voxel*, so in edema it fits an edema kernel and adapts — Step 6
prints how imperfectly, which is the honest version of that advantage. CSD is
handed a fixed response and cannot adapt. Setting `'matched'` hands CSD
something it could not have on real data.

**Measured, noise free, 60° crossing** (truth 60.00; band-limited truth 60.94):

| tissue | response used | SSST sep | MSMT sep |
|---|---|---|---|
| healthy | healthy | 60.94 | 60.94 |
| edema | **healthy** (the default) | 60.94 | **60.94** |
| edema | matched (edema) | 60.94 | **21.49** |

The result is the opposite of the intuition: giving MSMT the *matched* edema
response makes it much worse. The reason is conditioning, and it is measurable.
MSMT separates tissues on how the `l = 0` response decays across shells, and the
edema WM response decays almost like CSF:

| `l = 0` decay across shells | b=0 | b=1 | b=2 | b=3 |
|---|---|---|---|---|
| healthy WM | 1 | 0.513 | 0.335 | 0.253 |
| edema WM | 1 | 0.188 | 0.069 | 0.040 |
| CSF | 1 | 0.050 | 0.003 | 0.000 |

The angle between the WM and CSF response vectors falls from **31.3°** (healthy)
to **8.9°** (edema). At 8.9° those two columns of MSMT's design matrix are
nearly parallel, the WM/CSF split becomes unstable, and the WM fODF is
corrupted. So `'healthy'` is both the realistic choice and the better
conditioned one.

**One number to carry forward, because it is a tractography result rather than
an angular one.** With the healthy response on edema signal, MSMT assigns only
**4.7%** of the `l = 0` signal to WM (WM/GM/CSF = 0.047 / 0.285 / 0.668). The
peak *orientation* is perfect, but the fODF amplitude is tiny — which is exactly
the section 6.1 story: MSMT-CSD dims in edema. Good for terminating in CSF,
bad for tracking *through* an edematous region. Orientation accuracy alone will
not show this; the amplitude has to be reported with it.

### Where the response files are, and what sets the FOD order

`dwi2fod` needs three response `.txt` files for `msmt_csd` (WM, GM, CSF) and one
for `csd`. **They are generated at run time**, not tracked, into
`deconv_comparison/mrtrix/` — which is gitignored, because they are derived
entirely from the kernel and are reproduced exactly on the next run:

```
mrtrix/ms_<preset>_resp_wm_lmax<L>.txt      one row per shell, L/2+1 columns
mrtrix/ms_<preset>_resp_gm.txt              one row per shell, 1 column
mrtrix/ms_<preset>_resp_csf.txt             one row per shell, 1 column
mrtrix/ms_<preset>_resp_wm_b3_lmax<L>.txt   top shell only, for SSST-CSD
```

`<preset>` is `healthy` or `edema`, so both tissues keep their own set. They
persist after a run — read them with `shview` or `cat` to see exactly what each
arm was deconvolved with.

**`-lmax` is not passed to `dwi2fod`, and does not need to be.** MRtrix's
documented default is *"the lmax of the corresponding response function, based
on its number of coefficients, up to a maximum of 8"* — and one WM response is
written per Lmax with exactly `L/2+1` columns, so the file already pins the
order. The same mechanism gives `msmt_csd` its GM and CSF at lmax 0, since those
responses are one column wide.

Measured rather than assumed: with and without the flag, `dwi2fod` returns the
same **15, 28 and 45** coefficients at Lmax 4, 6 and 8, and the output FODs are
**bit identical** (`max|diff| = 0`).

**Images are written as `.mif`** — one self-contained file per image, rather than
the `.mih` + `.dat` pair the older pipeline scripts use. `mrtrix_io.m` picks the
format from the extension it is given, so the older callers are unchanged. The
writer resolves the `file: . <offset>` circularity (the offset is part of the
header whose length it describes) by iterating to a fixed point, and it is
verified against MRtrix on every run via `mrinfo`, plus a round trip through
`mrconvert` that matches the `.mih` path bit for bit.

### Findings, not failures

**RETRACTED: the large MSMT-CSD noise-free error was a setup error, not a
property of the method.** An earlier version of this file reported MSMT-CSD at
7.80° noise-free against a 1.27° ceiling and attributed it to response
mismatch, citing section 6.4. That was wrong, and the user was right to
disbelieve it — it contradicted Jeurissen et al. 2014.

**The cause.** `dwi2fod csd` and `dwi2fod msmt_csd` do not ship comparable
defaults. The SSST algorithm's non-negativity constraint has strength 1;
`msmt_csd`'s `-neg_lambda` defaults to **1e-10**, essentially unregularised.
Running both "at their defaults" compared a constrained arm against an
unconstrained one. The MSMT fODF came back much blunter — `l = 6` band power
**0.31** against SSST's **1.11** — which pulls a crossing's two lobes together
and displaces the peaks.

Peak separation on the noise-free crossing (true 60.00; band-limited truth at
Lmax 6 gives 60.94):

| true angle | SSST-CSD | MSMT, MRtrix defaults | MSMT, `-neg_lambda 1` |
|---|---|---|---|
| 60° | 60.94 | **48.67** | **60.94** |
| 75° | 73.95 | 70.99 | **73.95** |
| 90° | 89.36 | **86.38** | **89.36** |

MSMT under-separated at **every** angle including 90°, which no published
MSMT-CSD comparison shows — that was the tell that the setup was wrong rather
than the method. With `-neg_lambda 1` it matches SSST-CSD exactly at 60, 75 and
90 degrees.

**Excluded first, each by measurement rather than argument:** the shells
(`msmt_csd` on b = 3 alone fails identically), the tissue count (WM-only
identical), the response family (delta only partly helps), the peak finder
(`sh2peaks` agrees to 0.06°), the SH basis (SSST is exact), and signal scaling
(bit identical from S0 = 1 to 10⁴).

`MSMT_NEG_LAMBDA` and `MSMT_NORM_LAMBDA` are now explicit in the Configuration
block, and Step 6b prints the `l >= 2` band power of both arms on every run, so
a repeat of this failure is visible in the output rather than only in the peak
table. **Report both values with any MSMT number** — the sensitivity spans
48.67° to 64.97° at a true 60° crossing, so an MSMT result quoted without them
is not reproducible.

**Also disproved: section 6.4's low-b hypothesis.** The Monte Carlo report
supposed that low-b shells carry too little angular contrast and pull the joint
fit toward a single lobe. Measured here, SSST-CSD on **b = 1 alone** recovers
60.94° — identical to b = 3 alone. Low-b data is not the problem.

**The angular error uses the largest peak only, and a symmetric crossing is a
near-tie.** MSMT's two lobes differ by ~7% in amplitude; which is "primary" can
flip on the last few ulps, and the reported error then jumps between the two
lobes' errors with no real change in the reconstruction. Read it alongside the
peak count and the spurious count. The peak finder lives in
`helpers/fODF_peak_score.m` and records this in its header.

**Noise can improve a bias-dominated arm's mean angular error**, so the "more
noise must not help" check is classified rather than simply failed: with the
noise gone, what is left is whatever systematic displacement the method has, and
noise partly masks it by jittering the primary lobe around. The check reports an
arm whose noise-free error is more than twice its ceiling as a **NOTE**, and
fails only an arm that sits near its ceiling and still degrades without noise.

**The numbers that motivated this were the void MSMT ones** (see the retraction
above), so the classification has not yet been exercised by a case known to be
real. It is kept because it is right in principle and costs nothing, not because
it has been demonstrated.

**Sigma is recovered at b = 0, not over all volumes.** Rician noise is a
magnitude, so `std(S_noisy − S_clean)` equals sigma only where the signal is
well above the noise floor — and in edema it is not: the mean b = 3 signal is
**0.0399** against `sigma = 0.1`, i.e. two and a half times *below* the noise,
where the estimate reads 14% low. At b = 0 the signal is exactly 1 by
construction, so that is the one place the residual is clean Gaussian. Both
numbers are printed; the gap between them is the Rician floor.

### What the smoke run already shows

Indicative only — `NREP = 27` — but the ordering is the published one, and now
measured on **one noise draw shared by all three arms**. 60 degrees, Lmax 6:

| kernel | arm | SNR 10 correct / spurious | SNR 30 correct / spurious | SNR inf error |
|---|---|---|---|---|
| healthy | SMI | **100.0% / 0.000** | 100.0% / 0.000 | 1.72° |
| healthy | SSST-CSD | 44.4% / 0.667 | 100.0% / 0.000 | **1.27°** (the ceiling) |
| healthy | MSMT-CSD | 85.2% / 0.000 | 92.6% / 0.000 | 7.80° |
| edema | SMI | 0.0% / 1.741 | **70.4% / 0.296** | 1.27° |
| edema | SSST-CSD | 0.0% / 3.111 | 0.0% / 2.556 | 1.27° |
| edema | MSMT-CSD | 33.3% / 0.074 | 48.1% / 0.000 | 11.72° |

Two things stand out. **In healthy tissue the arms trade places with SNR** —
constrained SMI is the noise-robust one, SSST-CSD the high-SNR
angular-resolution one, which is exactly section 6.4 of the Monte Carlo report,
here reproduced without any cross-run comparison. **In edema the trade
disappears and SSST-CSD collapses**: 0% correct at both finite SNRs, with 2.5–3
spurious peaks per voxel, where SMI still recovers 70% at SNR 30. That is the
"less anisotropic signal, amplified in inverse proportion" mechanism showing up
as a method difference rather than as a claim.

## `smi_manuscript_60deg.ipynb` — the three-arm notebook

The same experiment as the `.m`, with the CSD arms wired in rather than stubbed,
and cut to the **healthy kernel only**; the edema arm waits on a synthetic
response.

| arm | what it is |
|---|---|
| SMI | `SMI.fit`, `fODF_regularization.flag_nonneg = 1`, `CS_phase = 0` |
| SSST-CSD | `dwi2fod csd` on the top shell — the MRtrix3 binary |
| MSMT-CSD | `dwi2fod msmt_csd` on all four shells, 3 tissues — the binary |

**All three responses are the exact analytic kernel response**,
`r_l(b) = K_l(b) sqrt((2l+1) 4 pi)`, so response estimation is not a confound
and a difference between arms is the deconvolution. Swapping in an estimated
response is one `RH.read_response` call in Step 6b — but see the warning there
about which acquisition a stored response was estimated on.

Three things the `.m` file's commented hook asserted that are **not true**, all
corrected in the notebook:

- **The arms do not "appear automatically" in the figures.** The scoring arrays
  and every figure were indexed by `(Lmax, SNR, condition)` with no arm
  dimension at all. There is one now.
- **The peak finder was not scale-free.** It subtracted the constant
  `1/(4*pi)`, which is the isotropic part only for an fODF with `p_00 = 1` in
  every voxel. An MRtrix FOD is unnormalised and its `l=0` term varies per
  voxel and carries apparent fibre density. The notebook subtracts **each
  voxel's own `l = 0` term**, and checks that this agrees with the constant for
  the SMI arm.
- **The response must be evaluated at MRtrix's own per-shell b**, read back with
  `mrinfo -shell_bvalues`, not at the nominal 0/1/2/3. This protocol's shells
  jitter — MRtrix reports `[0, 998.28, 1998.17, 2996.06]` s/mm² — so the
  nominal values attach each response row to a b nobody acquired at.

`SMOKE_TEST = true` is the shipped default: one Lmax, three SNRs, `NREP = 27`,
a few minutes. `false` gives the manuscript configuration.

**Its figures have been verified to render**, which is not true of anything else
in this repository. The Octave Jupyter kernel produces inline PNGs, so Figure 1's
shared radial scale is confirmed working rather than assumed. Only the `gnuplot`
toolkit is available in this container, and it ignores `camlight` and `lighting`,
so glyphs are flat-shaded here and correct in MATLAB. Each figure is wrapped in
`try`/`catch` and prints `** FIGURE FAILED **` rather than taking the run down.

`smi_manuscript_60deg.m` keeps **every knob in one Configuration block at the
top**, including the kernels, so retuning it never means opening another file.
It simulates **both tissues in one run**, back to back off the same ground
truth, protocol and noise seeds — only the kernel differs:

| preset | `[f Da Depar Deperp fw]` | extra-axonal |
|---|---|---|
| `healthy` | `[0.60 2.0 2.0 0.50 0.02]` | 0.38 |
| `edema` | `[0.10 2.4 2.7 1.15 0.35]` | 0.55 |

The `edema` kernel is shaped from a *fitted* kernel, which is why the
diffusivities are not round and why `Depar` exceeds `Da`. Two values depart from
the fit as supplied, both deliberately: **`f` is 0.10 rather than 0.05**, because
0.05 is exactly the lower bound of SMI's default training prior and a truth
sitting on the prior boundary is estimated with a one-sided bias that is
indistinguishable from a real effect; and **`fw = 0.35`**, which the fit did not
specify, modelling edema as added free water. `Deperp = 1.15` remains close to
its own prior cap of 1.2 — left as fitted, but the first place to look if the
extra-axonal diffusivities come back biased low.

Only three stateless utilities still come from `mc_config.m` — `pick_grid`,
`rotate_about`, `load_protocol_file`. They are shared rather than copied so the
fibre-axis convention and the protocol reader cannot drift from
`gen_montecarlo.m`. Changing a knob in the notebook does **not** change the
pipeline, which keeps its own values in `mc_config.m`.

`smi_simulation_walkthrough.m` draws four figures: the ground truth fibre
geometry, the kernel as a response function, the forward convolution with and
without noise, and the reconstructed fODFs at each Lmax.

`smi_manuscript_60deg.m` draws six, and every 3D panel opens in an **isometric**
view so all subfigures are readable without touching the camera:

| figure | layout |
|---|---|
| 1 | the ground truth fODF, then **each kernel's response glyph per shell**, then a **third row holding the difference between them**. All response glyphs share **one radial scale**, so size carries meaning: a smaller glyph is a genuinely smaller signal |
| 2 | the signal as surfaces on the sphere, **healthy**: b down the rows, SNR across |
| 3 | the same for **edema**, on the same scale as Figure 2 |
| 4 | **the SMI arm's** reconstruction, **healthy**: Lmax down the rows, truth in column 1 then one column per SNR |
| 5 | the same for **edema** |
| 6 | **bias, spread and spurious peak count against SNR** — one curve per Lmax, **one row per kernel AND arm**, y limits shared per column |
| 7 | **the three arms side by side**, at `LMAX_FIG`: one row per arm, truth then one column per SNR. One figure per kernel |

Figure 7 is the comparison figure. Each MRtrix FOD is put onto SMI's `p_00 = 1`
convention before the shared renderer sees it, by dividing out its own `l = 0`
term — a change of **scale only**, which multiplies every band by the same
factor and so preserves peak orientation exactly. Without it the SMI glyphs
would vanish next to the CSD ones, which live on an unnormalised amplitude
scale.

**Only Figure 1 uses a shared radial scale**, because it is the one panel whose
point is that the edema response is genuinely smaller than the healthy one.
Everywhere else each panel autoscales and fills its box: Figures 2 and 3 are
about the *shape* of the signal at each shell and SNR, Figures 4 and 5 about the
shape of the reconstruction, and shrinking one of them into illegibility would
buy nothing.

Getting that shared scale to work takes two things and only one is obvious.
Scaling a glyph's radius does nothing on its own, because `axis equal` fixes a
panel's aspect ratio but *not* its limits — MATLAB autoscales each subplot to
its own data, so a glyph half the size gets an axis range half as wide and is
drawn at exactly the same size. Figure 1 therefore scales every glyph to a
maximum radius of 1 **and** pins each panel to `GLYPH_LIM`. Its ground-truth
fODF is deliberately *not* on the response scale, since an fODF and a signal are
different quantities.

Set `MAKE_FIGURES = false` in that section to skip them. In MATLAB they appear
inline in the Live Script; under Octave graphics are often unavailable, and the
printed numbers are the deliverable either way — every figure is drawn from the
same arrays Step 7 prints, so the tables and the plots cannot disagree.

The commented hook that used to sit before Step 7 of `smi_manuscript_60deg.m`
is **gone**, replaced by Step 6b, which runs the arms rather than describing how
someone might. Three things it claimed were wrong, and all three had to be fixed
rather than carried over:

- **The arms do not "appear automatically" in the figures.** The scoring arrays
  and every figure were indexed by `(Lmax, SNR, condition)` with no arm
  dimension at all. There is one now, and it is the leading index.
- **The peak finder was not scale-free.** It subtracted the constant
  `1/(4*pi)`, which is the isotropic part only for an fODF with `p_00 = 1` in
  every voxel. An MRtrix FOD is unnormalised and its `l = 0` term varies per
  voxel. Step 7 now subtracts **each voxel's own `l = 0` term** and checks that
  this still equals the constant for the SMI arm.
- **The response must be evaluated at MRtrix's own per-shell b**, read back with
  `mrinfo -shell_bvalues`, not the nominal 0/1/2/3. This protocol's shells
  jitter — MRtrix reports `[0, 998.28, 1998.17, 2996.06]` s/mm² — so the
  nominal values would attach each response row to a b nobody acquired at.

The hook's fourth claim, that `dwi2fod` must be given `-lmax` matching
`LMAX_LIST` or the ceiling is the wrong bound, is **half right and no longer
how it works**: the order does have to match, but the response file's column
count already enforces it, so the flag is redundant. See above.

`SMOKE_TEST = true` cuts the file to one Lmax, three SNRs and `NREP = 27`, which
runs in minutes and still executes every CHECK. Its numbers are indicative only
and every printout says so. `false` is the shipped default and the manuscript
configuration: 42 `SMI.fit` calls, hours. The two MRtrix arms cost seconds
either way, because they scale with voxel count rather than with fit count.

### Testing the CSD arms on their own

`smi_manuscript_60deg.m` runs all three arms, but its SMI arm costs 42 fits and
hours at the manuscript settings, which makes it the wrong instrument for "does
the MRtrix side still work". From `deconv_comparison/`:

```
octave-cli --no-gui -q test_csd_arms.m          # ~2 seconds
```

It builds the same noise-free signal from the same kernel and protocol, hands it
to `dwi2fod`, and scores the peaks at 60, 75 and 90 degrees — the whole MRtrix
path, none of the cost. It asserts that SSST-CSD recovers the true angle to
within the direction grid, that MSMT-CSD separates as well as SSST-CSD, and that
both handle 90 degrees. That last one is the check that would have caught the
`-neg_lambda` bug immediately.

`test_csd_arms(true)` adds the sensitivity sweep:

| `neg_lambda` | `norm_lambda` | separation at a true 60° crossing |
|---|---|---|
| 1e-10 | 1e-10 | 48.67 ← **MRtrix's shipped defaults** |
| 1e-3 | 1e-10 | 53.65 |
| 1 | 1e-10 | 64.97 |
| 1e-3 | 1e-3 | 48.67 |
| **1** | **1e-3** | **60.94** ← what the manuscript file uses |

`neg_lambda = 1` is the principled half of that choice: it matches the strength
`dwi2fod csd` uses for its own non-negativity constraint, which is what makes
the two arms like-for-like. **`norm_lambda = 1e-3` is the least justified
constant in this package** — it was chosen because it lands on the band-limited
truth, which is uncomfortably close to fitting the answer. Interrogate it before
publishing anything that depends on it.

### Check it before you run it

`check_manuscript_static.m` answers "will this file reach the end" in seconds,
without simulating anything:

```
octave-cli --no-gui -q check_manuscript_static.m
```

It parses the whole file with `if false ... end` so every line is checked and
none executed, audits the **subscript arity** of the scoring arrays, and
verifies every `RUN{}` field that is read is also written.

The arity audit exists because of a specific bug. When the CSD arms went in,
`res_all` and friends gained a leading arm index and became 4-D; one read 500
lines below the change was missed, and a full run died in the final summary
table after everything else had already printed. The audit is **bracket-depth
aware** and **follows cell aliases** — the missed read was
`sums{im}(iL, SNR_ORD(k), ic)`, through a cell built from those arrays, so a
checker looking only for `res_all(` sails straight past it. That blind spot was
found by injecting the bug back into a copy and confirming the checker caught
it, which is the only way to know a checker works.

It does not tell you the numbers are right. It tells you the run will finish.

### Why the edema arm looks noisier, when the noise is identical

`sigma = 1/SNR` against `S0 = 1` in both runs, so the noise added is exactly the
same. What shrinks is the **anisotropic signal**. All three edema parameters push
the same way: `f` 0.60 → 0.10 removes most of the stick compartment, which is the
anisotropic one; `Deperp` 0.50 → 1.15 makes the extra-axonal tensor nearly
isotropic; and `fw` 0.02 → 0.35 adds purely isotropic signal. `K_2`, `K_4`, `K_6`
all shrink while `K_0` stays comparable — and the orientation information lives
entirely in those `l >= 2` invariants.

Deconvolution then amplifies rather than rescues, because recovering `p_lm` means
**dividing by `K_l`**: the gain is `1/K_l`, so smaller invariants amplify the same
noise more. Same mechanism as the CSF blow-ups documented elsewhere in this repo
(`g_2 = 1/||K_2||`, 3.6 in CSF against 0.6 in white matter).

Figure 1's shared scale is what makes this visible, and the table it prints gives
the numbers: on that scale the healthy response peaks at 0.83 / 0.74 / 0.68 for
b = 1 / 2 / 3, the edema response at 0.29 / 0.16 / 0.12, and both are 1.000 at
b = 0 by construction. So the framing is not
"edema is noisier" but **"edema has less anisotropic signal, and the
deconvolution amplifies noise in inverse proportion to it."**

## Running it

```
cd deconv_comparison/notebooks
octave-cli --no-gui -q smi_simulation_walkthrough.m      # about 6 minutes
```

or from MATLAB, from anywhere:

```matlab
run('deconv_comparison/notebooks/smi_simulation_walkthrough.m')
```

Nothing needs to be generated first and no data is downloaded. The acquisition
is a **real HCP 3-shell scheme**, read from `../protocol/hcp_real_3shell.txt`
via `mc_config.m`, so the MATLAB side of the package needs no Python.

Almost all the runtime is Step 6. In `smi_simulation_walkthrough.m` that is one
`SMI.fit` per Lmax, and `SMI.fit` is dominated by training its regression rather
than by voxel count, so each of the three Lmax values costs the same ~2 minutes
whether there are 100 voxels or 10,000.

**`smi_manuscript_60deg.m` is a different order of magnitude**, and at its
defaults it runs in hours rather than minutes. It fits
`NKERN × numel(SNR_LIST) × numel(LMAX_LIST)` times — **42** at two kernels,
`SNR_LIST = [5 10 20 30 50 100 Inf]` and `LMAX_LIST = [4 6 8]` — each on `NREP`
voxels, and at
`NREP = 1000` the per-voxel constrained deconvolution is no longer negligible
next to the training. **`NREP` is the knob**: drop it to 25 for a version that
runs while you read it. The file prints the elapsed time of every fit as it
goes, so the first two lines of Step 6 tell you what the whole run will cost.

One entry of the sweep is cheaper than it looks and dearer than it should be:
at `SNR = Inf` there is no noise, so all `NREP` realisations of a condition are
the same signal and come back bit identical. That block costs a full fit to
produce one distinct answer. It earns its place as the reference column of
Figures 3–5 and as the check that the fit reaches the truth at all, but it is
the cheapest entry to drop if runtime is what is hurting.

MRtrix is **not** required. The last step writes the fODFs as MRtrix SH images
so you can check them with `sh2peaks` or `mrview` yourself; the walkthrough
prints those commands and runs none of them.

## Opening it as a MATLAB Live Script

These are plain `.m` files written with publish-style markup (`%%` sections,
`*bold*`, `|monospace|`), so MATLAB converts them with no edits:

- **In the editor**: open the `.m`, then *Save As* and choose `.mlx`.
- **From the command line**:
  ```matlab
  matlab.internal.liveeditor.openAndSave( ...
      'smi_simulation_walkthrough.m', 'smi_simulation_walkthrough.mlx');
  ```

The `.m` is kept as the source of truth rather than the `.mlx` because `.mlx`
is zipped XML: it cannot be diffed, reviewed in a pull request, or run under
Octave, and this package's whole claim is that an outsider can check it.

## What "verifiable" means here

Every step ends with one or more `CHECK` lines that compare its output against
something computed a *different* way, and print `ok` or `** FAILED **` next to
the number. Reading only the `CHECK` lines is a complete audit. The ones worth
understanding:

| check | why it is not circular |
|---|---|
| zonal response == SMI forward model | the same single-fibre signal built from the response-function form and from `SMI`'s own convolution. This identity is what lets an SMI kernel be handed to MRtrix at all |
| harmonics vs direct convolution | the signal computed with no spherical harmonics anywhere, by summing over the quadrature grid. Catches a basis or normalisation error that a self-consistent SH round-trip would hide |
| spherical mean is orientation-blind | SMI fits the kernel from rotational invariants, so the `l = 0` invariant must not depend on where the fibres point |
| recovered sigma | the noise level read back out of the simulated data, rather than the value that was typed in. In the manuscript file this runs once per SNR |
| the ceiling | the same peak finder run on the noise-free truth *truncated to the same Lmax as the fit*. Separates "the method failed" from "this angular order cannot represent the answer" |
| SMI bins the shells | the real protocol has 18 distinct b values; `SMI.Group_dwi_in_shells_b_beta_TE` must recover `[18 90 90 90]`, or every rotational invariant downstream is built from a fraction of the directions |
| the noise-free arm reaches the truth | *(manuscript file only)* with `sigma = 0` the only things left between fit and truth are the band limit and the estimated kernel, so the `SNR = Inf` row must recover the true fibre count wherever the truth itself resolves. If it does not, no finite-SNR row below it is interpretable |
| the noise-free arm is deterministic | *(manuscript files only)* at `sigma = 0` every realisation of a condition is fed the same signal, so the block must come back the same. **See the correction below — the `.m` file states this invariant in a form that is false, and the `.ipynb` states the measured one** |
| more noise does not help | *(manuscript file only)* the mean angular error at the top of the sweep must not exceed the error at the bottom. Compared at the extremes rather than pairwise, so Monte Carlo error between adjacent SNRs cannot trip it |

In the manuscript file the ceiling is not scored by a *copy* of the peak finder:
the band-limited truth is prepended as the first row of every scored block, so
it goes through the same lines as the realisations. That is why the `SNR = Inf`
row and the `ceiling` row of each table print the same angular error.

### Correction: `SMI.fit` is not bit-reproducible voxel to voxel

`smi_manuscript_60deg.m:1069` checks that the noise-free block has standard
deviation **exactly** `0`, and the table above used to describe that as the
invariant. **Both are wrong**, in two separate ways, and the notebook's smoke
run found them. The `.m` file has never been run at full size, so neither had
ever been exercised.

1. **`std(x) == 0` is not a determinism test.** `std` computes `sum(x)/n` first,
   and for *n* identical values that division need not return the value exactly.
   A perfectly deterministic block can yield `std ~ 1e-16`. Whether it does
   depends on the bit pattern of the particular angular error — which is why
   this check failed for two arms and passed for a third on identical logic.

2. **`SMI.fit` genuinely is not bit-reproducible.** Measured directly: feed 27
   voxels on a `[3 3 3]` grid a bit-identical signal at `sigma = 0`, with the
   non-negativity constraint converging in the same 4 iterations in every
   voxel, and SMI returns **three distinct answers — one per slice** — differing
   by `1.6e-12` in `plm` and `1.1e-11` in the kernel. That is BLAS blocking:
   reduction order inside the matrix operations depends on where a voxel sits
   in the array. It is twelve orders of magnitude below the peak finder's 2.6°
   grid resolution and changes no result.

**Both MRtrix arms *are* bit identical**, because `dwi2fod` works one voxel at a
time. So the corrected check is worth keeping rather than deleting: it separates
a method that is exactly reproducible voxel to voxel from one that is
reproducible only to floating point. The notebook tests the SH coefficients
directly against a `1e-9` tolerance and prints the measured number.

Three results in the output are *findings, not failures*, and each is called out
in place because it looks alarming otherwise:

- **The band-limited ground truth is negative over roughly 40% of the sphere.**
  Truncating a Watson mixture at `Lmax` rings, and the rings cross zero. The
  non-negativity constraint is therefore a regularizer, not a statement of
  fact — the truth does not satisfy it either.
- **The 30 degree crossing is resolved 0% of the time, and that is correct.**
  The ceiling shows the noise-free truth also comes back as a single peak: two
  Watson populations at `kappa = 16` merge into one lobe well before any
  deconvolution is involved. This is a property of the *ground truth*, not of
  any method, and it is worth knowing where the boundary sits — measured on the
  truth at Lmax 8, a 30 degree crossing needs `kappa >= 48` (about 8 degrees of
  dispersion) before it separates at all, while 45 and 60 degrees separate at
  every `kappa` from 8 upward. `mc_config.m` keeps `kappa = 16` because real
  white matter disperses; the cost is that the 30 degree condition currently
  measures the band limit rather than the method.
- **A loud `WARNING` about gradient directions on every run.** The supplied
  `.bvec` is unit only to `1.1e-6`, because it was written at seven significant
  figures. `mc_config.m` warns and normalises. It is deliberately noisy: left
  uncorrected it degrades the zonal-response identity at Lmax 8 from `1e-15` to
  `5e-7`, and that failure is invisible unless something checks for it.

## The SNR sweep, and the one trap in it

`smi_manuscript_60deg.m` runs `SNR_LIST = [5 10 20 30 50 Inf]` as six
independent blocks of voxels, one `SMI.fit` per block per Lmax. Each block is
seeded separately (`mc_config.m`'s seed plus the SNR index) so a noise level is
reproducible on its own and adding or removing an entry does not silently
change the realisations of the others. The cost of that choice is that the
blocks do *not* share common random numbers, so a difference between two SNRs
carries the Monte Carlo error of both.

Four numbers are reported per cell, because a sweep is where they stop agreeing:
**correct fibre count** (what tractography consumes), **mean angular error**,
**its standard deviation**, and **spurious peaks per voxel**. A method can hold
its angular error while it starts inventing peaks, or lose the peak entirely and
report a confident wrong direction; no single column catches both.

**"Bias" here is the mean angular error, which is the mean of a non-negative
quantity.** It does not go to zero for a perfect estimator — the floor is set by
the direction grid (about 1.5° at 1500 directions) and by the band limit. Read
it against the noise-free point of the same curve, not against zero. That is
what the `SNR = Inf` column is for.

**The trap: SMI does not fit the kernel at the sigma you pass in.** It
normalises sigma by the measured `b = 0` signal, bins the result into `Nlevels`
equal bins spanning `sigma_norm_limits`, and trains one polynomial regression
per occupied bin, evaluated at the **bin centre** (`SMI.m:2222-2302`). With the
shipped `sigma_norm_limits = [0 0.2]` (`SMI.m:562`) and `Nlevels = 10`
(`SMI.m:388`) that has two consequences only a sweep makes visible:

- **two SNRs can land in the same bin** and then share a regression trained at
  one noise level — σ = 0.02 (SNR 50) and σ = 0.033 (SNR 30) both fall in bin 2;
- **an SNR below 5 is off the end of the trained range.** σ = 0.2 is exactly the
  top edge, so anything noisier is clamped into the top bin and fitted with a
  regression trained at a *lower* noise level than the data actually has.

Step 5 prints the bin each SNR will land in and flags any that are clamped. The
bins are nominal: sigma is normalised by the *measured* `b = 0` signal, so
voxels of one SNR can straddle an edge.

## Why the manuscript uses 60 degrees only

30 and 45 stay in the general walkthrough, but neither makes a clean figure.
Measured on the **noise-free ground truth**, with no fitting and no noise:

| angle | Lmax 4 | Lmax 6 | Lmax 8 | true |
|---|---|---|---|---|
| 30° | 1 | 1 | 1 | 2 |
| 45° | 1 | 2 | 2 | 2 |
| 60° | 2 | 2 | 2 | 2 |

At `kappa = 16` a 30 degree crossing does not separate at any angular order —
it needs `kappa >= 48`, about 8 degrees of dispersion. 45 separates only from
Lmax 6 upward. **60 separates everywhere**, which is what makes it the
configuration where a difference between methods is attributable to the method
rather than to the band limit.

## Relationship to the pipeline

The walkthrough is not a copy of `gen_montecarlo.m`; it reads its constants
from `../mc_config.m`, the same file the pipeline reads, so the two cannot
drift apart about what is being simulated. The SNR sweep keeps that
correspondence: `gen_montecarlo.m` takes one SNR per invocation, and the
manuscript file's Step 6 makes one `SMI.fit` call per SNR for the same reason —
no fit ever sees two noise levels at once. It does differ from the pipeline in
two deliberate ways: peaks are found on a fixed direction grid so the operation
is visible in the file, where the pipeline uses `sh2peaks` for every arm; and it
fits at three angular orders where the pipeline fixes Lmax 6 to keep the three
methods comparable. Step 8 exports the fODFs so `sh2peaks` can be run on them
directly.

The manuscript file's export covers the **whole sweep** in one image per Lmax —
every SNR block, in the order Step 5 laid them out — so a single `sh2peaks` run
reaches all noise levels. `export/voxel_key.txt` therefore carries an **SNR
column** alongside the crossing angle and the true fibre axes, and its header
records how many contiguous blocks of how many voxels are in the file, and in
what order. The noise-free arm appears there as `Inf`.

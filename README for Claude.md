# README for Claude

Handoff for the next agent working on SMI fODF tractography through edema.

This is the **third** version of this file. Each version corrects conclusions the
previous one stated as settled, so treat section 2 as the most important thing
here: it is the list of things a confident-sounding earlier handoff got wrong,
including things *this* session's agent got wrong and then measured properly.
Earlier versions are at `409e3ca` and `900a8b6~1` in git history.

**Read sections 1 and 2 before proposing anything.**

---

## 1. Status

The framework is built. The most recent session moved the comparison work off
reimplementations and onto MRtrix3 itself, and changed one shipped default.

### On master

| feature | option | default | documentation |
|---|---|---|---|
| Regularized deconvolution | `options.fODF_regularization` | off | `Reports/REPORT_fODF_regularization_sweep.md` |
| Anisotropy modulation *(semi-retired)* | `options.fODF_modulation` | off | `Archive/README.md` |
| Post hoc outlier cap | `options.fODF_outlier` | off | `Reports/REPORT_fODF_outlier_cap.md` |

All three are opt-in, and `out.plm`, `out.pl`, `out.kernel` are bit identical
whether any of them is on or off. That invariant is tested, not asserted — keep
it. Patches `Patches/0001`-`Patches/0016` record each change and apply with
`git am --3way`.

The shview-style response viewer, the `deconv_comparison/` Monte Carlo package,
`Reports/REPORT_SMI_deconvolution_MonteCarlo.md`, its two figures, and the
`lambda_nonneg` default change from 10 to 1 were merged in PR #7. The shared
Monte Carlo configuration was merged in PR #9. Section 6.4 justifies the
default change.

**Current pipeline posture:** keep anisotropy modulation available for
reproducibility, but do not present it as an active pipeline recommendation.
The zonal-harmonics response viewer is likewise a learning and convention-check
exercise. Both are indexed in `Archive/README.md`.

### The manuscript simulation, `deconv_comparison/notebooks/`

The active work. Two narrated notebooks, both plain `.m` with publish-style
markup so they run in MATLAB and Octave and convert to a Live Script unedited:

| notebook | what it is |
|---|---|
| `smi_simulation_walkthrough.m` | the general SMI arm. Single fibre plus 30/45/60 degrees, one SNR, Lmax 4/6/8, four figures |
| `smi_manuscript_60deg.m` | **the manuscript figure source, and now the whole comparison.** 60 degrees only, two kernels back to back, SNR swept, and **all three arms — SMI, SSST-CSD, MSMT-CSD — running inside this one file** |
| `smi_manuscript_60deg.ipynb` | an earlier three-arm version on the Octave Jupyter kernel (PR #22). Healthy kernel only, exact delta response. **Superseded by the `.m`**; kept because its figures are verified to render inline |

Both put every knob in one Configuration block at the top and end every step
with a `CHECK` that compares its output against something computed a different
way. `notebooks/README.md` lists what each check establishes and why it is not
circular; read that before changing either file.

`smi_manuscript_60deg.m` at its defaults is **42 `SMI.fit` calls** (2 kernels x
7 SNR x 3 Lmax) and runs in hours. **`SMOKE_TEST = true` is the knob** — one
Lmax, three SNRs, `NREP = 27`, minutes, every CHECK still executed.

### All three arms now run inside `smi_manuscript_60deg.m`

**This replaced a two-language arrangement, deliberately.** An earlier version
of this work put the CSD arms in a Python notebook and exported results back to
the `.m` file. It worked and the identity between the two forward models was
measured (noise-free signal to 6.7e-15), but the user judged that two kernels
would invite methodological questions, and they were right: it is now one file.

Step 6b calls `dwi2fod csd` and `dwi2fod msmt_csd` on **the same `S_noisy`
array `SMI.fit` was just given**. One simulation, one noise draw, one peak
finder, one scope. There is nothing to reconcile, which is the whole point.

Three things about it that matter more than the code:

1. **The response is dispersion matched** (`RESPONSE_MODE = 'dispersed'`):
   `r_l(b) = K_l(b) p_l^Watson(kappa) sqrt((2l+1) 4 pi)`. A fibre here is a
   Watson at `kappa = 16`, not a delta, so this is the response a perfect
   `dwi2response` would recover. At b = 3 it lands between the delta and the
   three estimators and accounts for **~69% of the `l = 2` gap** — most of what
   makes an estimated response blunt really is dispersion, as section 6.3
   supposed. `'delta'` switches back. Because the response is built from the
   current iteration's kernel, **the CSD arms run for the edema tissue too**;
   the old "the edema arm waits on a synthetic response" blocker was an artefact
   of using a stored estimated response.
2. **`SMOKE_TEST` is new and is the knob to reach for first.** `true` gives one
   Lmax, three SNRs, `NREP = 27` — minutes, every CHECK still executed,
   indicative numbers only. `false` is the manuscript configuration: 42
   `SMI.fit` calls, hours. The MRtrix arms cost seconds either way.
3. **Step 7 gained a leading arm index** and its peak finder now subtracts each
   voxel's own `l = 0` term rather than the constant `1/(4*pi)` — mandatory,
   because an MRtrix FOD is unnormalised. The finder itself moved to
   `helpers/fODF_peak_score.m` so nothing can score an fODF a second way.

**Run `check_manuscript_static.m` before starting a long run.** It parses the
whole file without executing it, audits the subscript arity of the scoring
arrays, and checks every `RUN{}` field read is written — seconds, not hours.
It exists because adding the arm index made those arrays 4-D and one read 500
lines away was missed, killing a run in the final summary table after
everything else had printed. The audit is bracket-depth aware and **follows
cell aliases**, since the missed read was `sums{im}(iL, SNR_ORD(k), ic)`;
a checker looking only for `res_all(` misses it. That blind spot was found by
injecting the bug into a copy and confirming the checker caught it — do the
same to any check you add here, or you do not know it works.

**Two CHECKs had to be corrected, both found by running rather than reading:**

- *Sigma recovery moved to the b = 0 volumes.* Rician noise is a magnitude, so
  `std(S_noisy - S_clean)` equals sigma only where signal exceeds the noise
  floor. In edema the mean b = 3 signal is **0.0399** against `sigma = 0.1` —
  the signal is *below* the noise — and the old all-volume estimate read 14%
  low and failed. At b = 0 the signal is exactly 1 by construction. Both
  numbers print; the gap is the Rician floor.
- *"More noise must not help" is now classified, not just failed.* Its premise
  is that error is noise-dominated, which is false for an arm whose error is
  dominated by systematic bias: noise can jitter a displaced primary lobe closer
  to a true axis, so removing it makes the mean angular error worse. An arm
  whose noise-free error exceeds twice its ceiling is reported as a NOTE; only an
  arm sitting near its ceiling that still degrades without noise fails.
  **The numbers that originally motivated this were the void MSMT ones** (see
  the retraction below), so the check has not yet been exercised by a case known
  to be real. Keep it -- the classification is right in principle and costs
  nothing -- but do not cite it as evidence of anything.

**RETRACTED, and it is the most important correction in this section.** An
earlier version of this handoff reported MSMT-CSD's noise-free 60 degree error
as **7.80 deg against a 1.27 deg ceiling** and attributed it to response
mismatch, citing section 6.4. The user disbelieved it because it contradicts
Jeurissen et al. 2014, and they were right.

**It was a setup error.** `dwi2fod csd` and `dwi2fod msmt_csd` do not ship
comparable defaults: the SSST constraint has strength 1, while `msmt_csd`'s
`-neg_lambda` defaults to **1e-10**, essentially unregularised. Running both
"at their defaults" compared a constrained arm against an unconstrained one.
The MSMT fODF came back much blunter -- `l = 6` band power 0.31 against SSST's
1.11 -- pulling a crossing's lobes together.

| true angle | SSST-CSD | MSMT, MRtrix defaults | MSMT, `-neg_lambda 1` |
|---|---|---|---|
| 60 deg | 60.94 | **48.67** | **60.94** |
| 75 deg | 73.95 | 70.99 | **73.95** |
| 90 deg | 89.36 | **86.38** | **89.36** |

MSMT under-separated at *every* angle including 90 degrees. That is the tell:
no published MSMT-CSD comparison shows a 90 degree crossing being missed, so the
setup was wrong rather than the method. Ruled out first, each by measurement:
the shells (msmt on b = 3 alone fails identically), the tissue count (WM-only
identical), the response family, the peak finder (`sh2peaks` agrees to 0.06
deg), the SH basis (SSST is exact) and signal scaling (bit identical from
S0 = 1 to 1e4).

`MSMT_NEG_LAMBDA` and `MSMT_NORM_LAMBDA` are now explicit knobs, and Step 6b
prints both arms' `l >= 2` band power every run so a recurrence is visible in
the output. **Any MSMT number must be reported with both values** -- the
sensitivity spans 48.67 to 64.97 deg at a true 60 deg crossing.

**Section 6.4's low-b hypothesis is disproved.** It supposed low-b shells carry
too little angular contrast and pull the joint fit to a single lobe. Measured:
SSST-CSD on **b = 1 alone** recovers 60.94 deg, identical to b = 3 alone.

**The lesson, since this is the second time.** Both times a result "looked
wrong" to the user it was a mistuned regularisation default -- `lambda_nonneg`
last time (section 2.1), `-neg_lambda` this time. When two arms are compared
"at their defaults", check that the defaults are comparable before believing
the difference.

**The old commented hook before Step 7 is gone.** Three of its claims were
false: the arms do not appear automatically in the figures (there was no arm
dimension), the peak finder was not scale-free, and the response must be
evaluated at MRtrix's own per-shell b rather than the nominal 0/1/2/3.

### Still on a branch, not merged

**`claude/freewater-simulations` (`1be06ea`, 3 commits)** — the free-water
comparison package, `REPORT_fODF_freewater.md`,
`REPORT_fODF_compartment_split.md`, two figures. **Still no PR, still unmerged,
now three sessions old.** Every measured number quoted in
`Reports/REPORT_fODF_outlier_cap.md` sections 4-5 comes from this branch, so it should
land or those numbers lose their provenance. `deconv_comparison/binio.{m,py}`
and `kernel.py` are *copies* of files on that branch — de-duplicate if it lands.

### Never run on real data

Everything. No result in any report has touched a patient scan — including the
"real data" that the Monte Carlo response functions are estimated from, which is
a synthetic phantom (`deconv_comparison/gen_phantom.m`). Two driver scripts live
only on the user's machine and arrive as file uploads: `run_smi_hcp.m` (single
subject) and `run_smi_batch_mod.m` (5-subject edema batch). `run_smi_batch_mod.m`
runs `Lmax = [0 2 2 2 4 4]`, i.e. **max(Lmax) = 4**, which changes the ceiling
arithmetic — see section 3, item 5.

---

## 2. Corrections to earlier handoffs

### 2.1 Now resolved: `lambda_nonneg`

The previous handoff's item 3 said the `lambda_nonneg = 10` default "wants a
cross-check" and was "the single most concrete open task in the repo". **That is
done. The default is now 1** (`SMI.m:978`).

The two sweeps genuinely disagree, and the reason is that they score different
things:

- `examples/example_fODF_regularization_sweep.m` minimises **relative L2 error over the
  whole sphere**. That integral is dominated by the isotropic part and by
  negative mass, both of which a heavier penalty fixes. It prefers 10.
- The Monte Carlo scores **peak orientation, fibre count and `l >= 2` ACC**.
  Those are what tractography consumes. It prefers 1, at every SNR and every
  crossing angle.

Nothing above 1 buys further spurious-peak suppression — every constrained
setting sits at exactly 0.000 spurious peaks — and everything above 1 costs
angular resolution. The doc comment at `SMI.m:935` states both measurements so
the next reader does not rediscover the disagreement.

### 2.2 Corrections to *this* session's own earlier claims

Both were made against dipy and both were wrong. MRtrix settled them.

1. **The `CS_phase = 1` discrepancy is a 180 degree rotation about z, not a
   mirror through `z = 0`.** x and y are negated, z is unchanged. Measured with
   MRtrix's own `sh2amp` + `sh2peaks`: peak at `[-0.565, 0.942, 1.527]` against a
   truth of `[0.565, -0.942, 1.527]`, 71.50 deg error.
2. **There is no extra `sqrt(2)` at `m != 0`.** That was dipy's `tournier07`
   basis differing from MRtrix's, not SMI differing from MRtrix. At
   `CS_phase = 0`, **SMI's basis is MRtrix's exactly** — the map is the
   identity and the peak is recovered to machine precision.

`Reports/REPORT_SMI_deconvolution_MonteCarlo.md` carries the corrected version.
`deconv_comparison/check_mrtrix_basis.sh` re-measures it in about a minute
against the real binaries; run that rather than trusting either claim.

### 2.3 Still standing from the previous handoff

1. **`lambda_tikhonov` does not damage the high `l` bands, and does not do much
   of anything.**[^regularization-review-posture] Re-confirmed this session: 0.3 vs 0 moves the 45 degree error
   from 21.28 to 21.27 deg. `Reports/REPORT_fODF_modulation.md` §3 attributes the `pl4`
   loss to it and is **still wrong** and still unfixed. The real cause is the
   non-negativity constraint plus error in the estimated kernel, whose `K_l` at
   high `l` is tiny and very sensitive.
2. **The regularizer that is load-bearing against blow-ups is non-negativity,
   not Tikhonov.** With the constraint off at SNR 15 the CSF peak reaches 1.87
   and white matter is *dimmer* than CSF.
3. **§7.5's rejection of tissue-fraction weights is bounded, not universal.**
   A tissue-fraction weight fails on axonal loss and is safe under water
   redistribution. Which regime real peritumoral tissue is in decides whether
   `f` is usable, and that is measurable in a known ROI. It was historically
   the most promising untested weight, but modulation is now semi-retired; do
   not prioritize it unless a concrete real-data question reactivates the work.

[^regularization-review-posture]: **Reviewer posture:** We are increasingly
    skeptical of departing from established CSD/deconvolution regularization
    conventions, including nonstandard uses or interpretations of Tikhonov
    damping. A reviewer is likely to ask why a deviation is necessary. Begin
    from the conventional baseline and require a reviewer-facing motivation,
    an ablation, and evidence that the change improves the downstream quantity
    that matters.

---

## 3. Conventions that must not be got wrong

Pinned analytically and verified numerically (round-trip `max|err| = 1.2e-15`,
`tests/test_SMI_response_helpers.m`).

`out.plm` is in the **normalized** convention `p_00 = 1`, covering `l = 2..Lmax`
only (the `l=0` term is not stored).

```
SH coefficients :  f_lm = plm .* sqrt((2l+1)/(4*pi))
l = 0 term      :  f_00 = 1/sqrt(4*pi)        (constant in every voxel)
fODF amplitude  :  A(u) = 1/(4*pi) + sum_{l>=2,m} f_lm Y_lm(u)
Lmax from count :  Lmax = sqrt(2*Nlm + 9/4) - 3/2
```

Consequences that keep mattering:

1. **The fODF integrates to 1 in every voxel.** It carries no density
   information at all. A CSF voxel and a coherent WM voxel have equal mass.
2. **The isotropic floor is a fixed `1/(4*pi) = 0.0796`**, above MRtrix's
   default `iFOD2 -cutoff 0.05`. An unmodulated SMI fODF therefore passes the
   tractography termination test **in every voxel of the brain**, CSF included.
3. To rescale a peak to target `T`, scale the anisotropic part by
   `s = (T - 1/(4pi))/(peak - 1/(4pi))`, never `T/peak`. Both the modulation and
   the outlier cap depend on this.
4. Scaling all `l>=2` coefficients uniformly preserves peak orientation exactly.
   Clipping SH coefficients **independently** does not and can swing peak
   orientations. Never do the latter.
5. **The hard physical maximum of a band-limited fODF is
   `sum_{l even <= Lmax} (2l+1)/(4*pi)`**: 1.194 at Lmax 4, 2.228 at Lmax 6,
   3.581 at Lmax 8. Any amplitude ceiling should be compared against this and
   against what SMI actually recovers, which is much lower.

**The kernel *is* a response function, and the conversion is one line.** SMI's
forward model (`SMI.m:818-820`) is
`S(u)/S0 = sum_lm K_l(b) p_lm Y_lm(u) sqrt((2l+1)*4*pi)`, so the zonal harmonic
response — exactly the numbers an MRtrix response `.txt` file holds, one row per
shell — is

```
r_l(b) = K_l(b) * sqrt((2l+1) * 4*pi)
```

`helpers/SMI_response_helpers.m` implements this and `examples/example_SMI_response_shview.m`
checks it against SMI's own forward model on every run, aborting above 1e-10
(measured 1.6e-15). Use those rather than rederiving.

**`CS_phase` decides whether SMI's SH basis matches MRtrix's.** At
`CS_phase = 0` they are identical. At the SMI default of 1 they differ by
`diag((-1)^m)`, which rotates every fODF 180 degrees about z. Anything that
hands SMI coefficients to an MRtrix tool — `sh2peaks`, `tckgen`, `shview` —
must fit at `CS_phase = 0` or apply the sign flip. See section 2.2.

**A shared radial scale on a glyph figure needs TWO things.** This cost a whole
review cycle. `axis equal` fixes a panel's aspect *ratio*, not its *limits* —
MATLAB autoscales each subplot to its own data, so a glyph drawn at half the
radius gets an axis range half as wide and lands on the page at exactly the same
size. Scaling the radius alone is invisible. Scale the glyph **and** pin the
panel (`set(gca,'XLim',...,'YLim',...,'ZLim',...)`).

Corollary worth keeping: only use a shared scale where size is the point. In
`smi_manuscript_60deg.m` that is Figure 1 alone, where the whole claim is that
the edema response is genuinely smaller. Everywhere else a shared scale just
shrinks one panel into illegibility, and the figures deliberately autoscale.

**Two different p2 exist and they are not the same quantity.** Same for p4.

- `out.kernel(:,:,:,ip2)` — fitted jointly with the kernel by polynomial
  regression on rotational invariants. `ip2 = 6` without T2 fitting, `8` with.
  `SMI.grab_kernel_pl` resolves the ambiguity from `out.shells` row 4 (TE).
- `out.pl(:,:,:,1)` — the l=2 invariant of the **deconvolved** fODF, a raw
  Euclidean norm, unclipped.

With `Lmax = [0 8 8 8]`, `RotInv_Lmax` defaults to 4, so **`out.kernel` has 7
maps, not 6**: `[f Da Depar Deperp fw p2 p4]`.

---

## 4. Environment and validation

Linux container, **GNU Octave 8/9, no MATLAB**. Working setup:

```
apt-get update && apt-get install -y octave octave-statistics octave-image
apt-get install -y mrtrix3          # 3.0.4; needed for deconv_comparison/
pkg load statistics; pkg load image;
```

**MRtrix3 3.0.4 is now a real dependency** of `deconv_comparison/`, not a
convenience. `dwi2response`, `dwi2fod`, `sh2peaks`, `sh2amp`, `dwiextract`,
`mrconvert`, `mrinfo` and `mrdump` are all called directly. Nothing in that
package reimplements an MRtrix algorithm any more (section 7).

`SMI.fit` **does run end-to-end in Octave** with three shims. They are now in
the repo at `deconv_comparison/stubs/` — `oct_path.m` puts them on the path:

| shim | why |
|---|---|
| `round(x,n)` | MATLAB two-arg form, used in `Group_dwi_in_shells_b_beta_TE` |
| `discretize(x,edges)` | not in Octave, used in `StandardModel_MLfit_RotInvs` |
| `datetime()` | not in Octave, used only by the log writer |

Other traps, all hit at least once:

- **Octave cannot call functions defined at the end of a script**, MATLAB
  requires them there. That is why `helpers/SMI_response_helpers.m` and
  `helpers/fODF_modulation_helpers.m` are separate files returning function-handle
  structs rather than local functions.
- **`SMI.vectorize` takes a different branch if any spatial dimension is a
  singleton.** Always build simulation volumes with all three dims > 1.
- **`end` inside a comma-list expansion** (`fod(p{:}, 2:end)` where `p` is a
  `num2cell` result) is miscomputed in Octave. Use explicit indices.
- **`median(...,4,'omitnan')`** is version-dependent. `SMI.neighbour_median`
  deliberately takes the median by sorting instead (NaN sorts last in both), so
  it needs no nanflag and no toolbox.
- Graphics are limited, but **no longer unusable — this changed and the old
  claim is wrong**. With `octave` (not just `octave-cli`) installed, the gnuplot
  toolkit renders headless and `print -dpng` writes real files:
  a figure script driving `print -dpng` produces real PNGs without error. gnuplot
  still ignores `camlight` and `lighting`, so glyphs are flat shaded here and
  correct in MATLAB. `deconv_comparison/octave_test_stubs/` stubs
  `figure/histogram/subplot/...` for the older plotting examples that only need
  to run headless; prefer actually rendering now that it works.
- `strel('disk',r,0)` — Octave requires the explicit `N`; `N=0` is an exact
  disk and is the better choice in MATLAB too.
- Python needs only **`numpy matplotlib`**, and only for the older
  `run_all.sh` campaign — `smi_manuscript_60deg.m` needs no Python at all.
  `setup_protocol.py` additionally needs `dipy` for its gradient-direction
  repulsion. **cvxpy is no longer needed** — it was only there for the removed
  dipy MSMT-CSD path.
- **If you do need Python here, check which `python3` you get.** In this
  container `/usr/local/bin/python3` is 3.11 while `apt` installs
  `python3-numpy` for **3.12**, so a plain `python3 -c "import numpy"` fails
  with the misleading *"you should not try to import numpy from its source
  directory"*. The real error underneath is
  `No module named 'numpy.core._multiarray_umath'`, i.e. an ABI mismatch.
  `/usr/bin/python3.12` is the working interpreter.
- **A real `.bvec` is not unit to machine precision.** The supplied HCP file is
  unit only to `1.1e-6`, being written at seven significant figures. Anything
  that reads `g(3)` as `cos(theta)` without normalising inherits that: measured
  at Lmax 8 it degrades the zonal-response identity from `1e-15` to `5e-7`.
  `mc_config.m`'s `load_protocol` warns loudly and normalises, and it is the one
  reader so nothing can bypass it.
- **A real `.bval` records EFFECTIVE b.** The "b = 0" volumes are `b = 5 s/mm^2`
  and carry meaningless unit direction vectors; the shells jitter (18 distinct
  values across 4 shells). `SMI.Group_dwi_in_shells_b_beta_TE` bins them
  correctly — the notebooks check it recovers `[18 90 90 90]` — and
  `B0_SNAP` in `smi_manuscript_60deg.m` sets the near-zero shell to exactly 0 so
  `S(0)/S0 = 1` holds exactly.
- **`pkill -f <pattern>` matches the invoking shell** and will kill your own
  Bash tool call (exit 144). Use `pkill -x`. Same failure mode makes
  `until ! pgrep -f script.py` an infinite loop; poll for output files instead.

**Parse-checking is not enough.** Runtime errors reached the user twice that
way, because the interesting code sits inside `try` blocks. Two patterns that
do work:

1. Extract the post-fit section of the driver script, stub the NIfTI I/O, feed a
   synthetic `out`, run all flag combinations.
2. For structural edits, a **bracket-depth-aware** keyword balance check.
   A naive one is useless: `Pp(2:end-1,...)` and `addpath(..., '-end')` both
   produce false positives.

**Do not commit `deconv_comparison/data/` or `deconv_comparison/mrtrix/`.**
Both are gitignored now. A 10 MB `.dat` slipped through once and had to be
stripped from branch history with `git filter-branch --index-filter` and
force-pushed.

---

## 5. What is in the repo, and what each thing is for

| file | what |
|---|---|
| `SMI.m` | the toolbox. All three opt-in features live here as static methods |
| `README.md` | user-facing documentation of the recommended pipeline, with short status notes for archived work |
| `Archive/README.md` | learning exercises and semi-retired exploratory work; modulation and the zonal-harmonics viewer are indexed here |
| `examples/example.m`, `examples/example_SMI_SSM.m` | the original data-fit and sensitivity-specificity examples |
| `Reports/REPORT_fODF_regularization_sweep.md` | the original `lambda_nonneg` measurement (now superseded on the default; see section 2.1) |
| `Reports/REPORT_fODF_modulation.md` | the anisotropy weight measurement. **§3 is wrong on Tikhonov, see section 2.3** |
| `Reports/REPORT_fODF_outlier_cap.md` | the cap measurement |
| `examples/example_fODF_regularization*.m` | regularization examples and the sweep |
| `examples/example_fODF_modulation.m` + `helpers/fODF_modulation_helpers.m` | semi-retired 7-class learning exercise; retained for reproducibility, not current pipeline guidance |
| `tests/test_fODF_outlier_cap.m`, `tests/test_SMI_outlier_cap.m` | the cap's tests. Self-contained |
| `Patches/0001`-`Patches/0010*.patch` | one patch per measured change, `git am --3way`-able |

The response/deconvolution work now on master:

| file | what |
|---|---|
| `helpers/SMI_response_helpers.m` | archived learning helper: kernel → zonal harmonics → glyph |
| `examples/example_SMI_response_shview.m` | archived learning exercise for response conventions and visualization |
| `tests/test_SMI_response_helpers.m` | 8 tests, all passing under Octave |
| `deconv_comparison/` | the Monte Carlo package. Has its own `README.md` — read that one for the run order |
| `deconv_comparison/notebooks/` | the two narrated notebooks, and a `README.md` listing every `CHECK` and what it establishes |
| `deconv_comparison/protocol/` | `hcp_real_3shell.txt` is the acquisition in use; `hcp_like_3shell.txt` is the synthetic scheme the published tables were measured on, marked superseded |
| `Reports/REPORT_SMI_deconvolution_MonteCarlo.md` | ~580 lines, 10 sections, plus "Findings, shortest form" at the top |
| `Reports/deconv_tables.md` | every result table, generated by `deconv_comparison/tables.py` |
| `Figures/fodf_response_shview.png`, `Figures/fodf_deconv_montecarlo.png` | the two figures |

The division of labour in `deconv_comparison/` is deliberate and is what the
user asked for: **Octave/SMI generates the signals and fits the SMI arm; MRtrix3
does CSD, MSMT-CSD and every peak extraction; Python does bookkeeping, scoring
and figures.** The only MRtrix behaviour implemented locally is reading and
writing its image format (`mrtrix_io.{m,py}`), verified against `mrinfo` and
`mrconvert -strides`.

On `claude/freewater-simulations` only: `freewater_comparison/`,
`REPORT_fODF_freewater.md`, `REPORT_fODF_compartment_split.md`, two figures.

---

## 6. What the four measurement campaigns found

### 6.1 SMI vs CSD vs MSMT-CSD in free water

Same synthesised DWI fed to all three, one shared evaluation sphere, one shared
peak finder, all at Lmax 6.

**SMI's fODF amplitude is blind to free water; CSD's and MSMT-CSD's scale with
the tissue fraction.** Peak ratio at `fw = 0.40` over `fw = 0`, dilution model:
SMI 0.94-1.01, CSD 0.60, MSMT 0.62. Structural, not a tuning artefact — SMI's
fODF has unit mass by construction.

The flip side, and the whole tractography problem in one line: SMI leaves CSF at
**0.32** (four times its isotropic floor, 100% above MRtrix's cutoff) where
MSMT-CSD leaves it at **0.028** and 0%. **SMI cannot dim in edema and cannot
terminate in CSF; MSMT-CSD does both.** On orientation SMI matches MSMT-CSD and
beats CSD.

### 6.2 Edema as redistribution, not dilution

Hold `f = 0.60` and convert extra-axonal water to free water. This matters more
than the free water fraction does: the same `fw = 0.40` costs CSD 13% and MSMT
23% of amplitude instead of 40%, because an AFD-style fODF tracks *fibre*
content and `f` is what did not move. See section 2.3, item 3.

Also: **MSMT-CSD's CSF volume fraction is not a usable free-water estimate under
this model** — 0.019 where the truth is 0.20. SMI's `fw` reads 0.151.

### 6.3 Response estimation

Edema with intact axons is **not distinguishable** from healthy WM on a high-b
anisotropy criterion (b=3 anisotropy 0.677 vs 0.682), so it is enriched in the
selection — but the estimated response is unchanged to three decimals either
way. The contamination is benign because the contaminant looks like the thing
being estimated.

What is **not** benign: any estimated response is far blunter than an idealised
delta response, because it absorbs fibre dispersion. Confirmed against real
MRtrix this session — normalised to `l = 0` at b = 3:

| response | `r_2` | `r_4` | `r_6` |
|---|---|---|---|
| `dwi2response dhollander` | -0.703 | 0.342 | -0.107 |
| `dwi2response tournier` | -0.714 | 0.349 | -0.116 |
| `dwi2response fa` | -0.718 | 0.352 | -0.121 |
| exact delta kernel | -0.829 | 0.442 | -0.179 |

Every estimator is 15-40% blunter than the truth, and the three agree with each
other far more closely than any of them agrees with the kernel.

### 6.4 The Monte Carlo comparison (this session)

10,000 realisations per condition per SNR. Conditions: single fibre and 15/45/60
degree crossings, Watson `kappa = 16`, kernel `[0.60 2.0 2.0 0.50 0.02]`,
Lmax 6 fitting against Lmax 8 truth, `CS_phase = 0`. **All peaks from
`sh2peaks`, including SMI's** — peak extraction is identical across arms.

**SNR 50, angular error and fraction resolved:**

| | single | 15 deg | 45 deg | 60 deg |
|---|---|---|---|---|
| SMI constrained | 0.27 deg | — | 81.1% / 6.52 deg | 1.43 deg |
| SSST-CSD | 0.46 deg | — | 96.9% / 2.84 deg | 1.61 deg |
| MSMT-CSD | 0.35 deg | — | **0.0%** / 22.25 deg | 6.85 deg |

**At low SNR the ordering inverts.** SNR 10, 60 degrees resolved: SMI 99.8%,
SSST 49.6%, MSMT 87.2%; spurious peaks SSST 0.633 against SMI 0.002. SNR 5,
single fibre correct count: SMI 96.8%, SSST 15.6% (1.602 spurious peaks),
MSMT 100%. **Constrained SMI is the noise-robust arm and SSST-CSD is the
high-SNR angular-resolution arm.**

**The MSMT-CSD 45 degree failure is real and is not a bug in the setup.** The
user's instinct that a 45 degree failure looked wrong was half right — most of
*SMI's* 45 degree error was the mistuned `lambda_nonneg = 10` default, which is
why it changed. MSMT's was controlled two ways and survived both: **0.0% with
the exact kernel response** and **0.0% at Lmax 8**. Its 60 degree error, by
contrast, *is* response-limited (6.85 → 2.36 deg with the exact response, 2.50
at Lmax 8). Hypothesis, **not measured**: the low-b shells carry almost no
`l >= 4` contrast at 45 degrees and pull the joint fit towards a single lobe.
Testing that is a shell-weighting experiment — see section 8.

**The `lambda_nonneg` sweep** (SNR 30, 2000 reps, `sh2peaks`), which is what
changed the default:

| | constraint off | λ=1 | λ=3 | λ=10 |
|---|---|---|---|---|
| 45 deg resolved | 95.9% | 55.4% | 0.2% | 0.0% |
| ACC, single | 0.9800 | **0.9800** | | 0.9298 |
| ACC, 45 deg | | **0.9656** | | 0.8869 |

λ=1 has the best ACC at every condition. Spurious peaks are 0.000 at every
constrained setting, so nothing above 1 buys anything.

---

### 6.5 What the manuscript simulation is set up to measure

Not yet run at full size. What is already established, all measured on the
*noise-free ground truth* with no fitting and no noise, so these are properties
of the truth rather than of any method:

**Which crossings are resolvable at all** (peaks found in the truth truncated to
each Lmax, `kappa = 16`):

| angle | Lmax 4 | Lmax 6 | Lmax 8 | true |
|---|---|---|---|---|
| 30 deg | 1 | 1 | 1 | 2 |
| 45 deg | 1 | 2 | 2 | 2 |
| 60 deg | 2 | 2 | 2 | 2 |

**A 30 degree crossing never separates at `kappa = 16`** — it needs
`kappa >= 48`, about 8 degrees of dispersion. 45 separates only from Lmax 6 up.
**60 separates everywhere**, which is why the manuscript uses it: a difference
between methods there is attributable to the method and not to the band limit.
Note the consequence for `run_smi_batch_mod.m`, which runs `max(Lmax) = 4` on
real subjects: at that order a 45 degree crossing cannot be resolved at any SNR.

**The two kernels**, and why the edema one is what it is:

| kernel | `[f Da Depar Deperp fw]` | extra-axonal |
|---|---|---|
| healthy | `[0.60 2.0 2.0 0.50 0.02]` | 0.38 |
| edema | `[0.10 2.4 2.7 1.15 0.35]` | 0.55 |

The edema kernel is shaped from a **fitted** kernel, which is why the
diffusivities are not round and why `Depar` exceeds `Da`. Two values depart from
the fit deliberately: `f` is **0.10 rather than 0.05**, because 0.05 is exactly
the lower bound of SMI's default training prior and a truth sitting on the prior
boundary is estimated with a one-sided bias indistinguishable from a real
effect; and `fw = 0.35`, which the fit did not specify. `Deperp = 1.15` is still
close to its own prior cap of 1.2 — left as fitted, and the first place to look
if the extra-axonal diffusivities come back biased low.

**Why the edema arm will look noisier when the noise is identical.** It is not
noisier: `sigma = 1/SNR` against `S0 = 1` in both. The *anisotropic signal* is
3-6x smaller, measured on Figure 1's shared scale — healthy peaks at
0.83 / 0.74 / 0.68 for b = 1 / 2 / 3, edema at 0.29 / 0.16 / 0.12, both 1.000 at
b = 0 by construction. `f` 0.60 -> 0.10 removes most of the stick compartment,
`Deperp` 0.50 -> 1.15 makes the extra-axonal tensor nearly isotropic, and
`fw` 0.02 -> 0.35 adds pure isotropy. `K_2`, `K_4`, `K_6` all shrink while `K_0`
does not, and deconvolution **divides by `K_l`** — the same `g_2 = 1/||K_2||`
mechanism as the CSF blow-ups in section 6.3. Say "less anisotropic signal,
amplified in inverse proportion", never "noisier".

---

## 7. Dead ends, with evidence

Do not re-walk these.

- **Do not reimplement MRtrix.** An earlier version of `deconv_comparison/`
  ran CSD and MSMT-CSD through dipy and reimplemented `dwi2response dhollander`,
  `mrthreshold` and `amp2response` from the MRtrix source — 10 files, all now
  deleted. It cost a day and produced two wrong conclusions (section 2.2) plus a
  sign error in the `amp2response` monotonicity constraint that silently
  collapsed the WM response to isotropic. `apt-get install mrtrix3` and call the
  binaries. Its *conclusions* about MSMT-CSD did reproduce under the real thing;
  `git log` has the reimplementation if it is ever wanted back.
- **`p2product` (the shipped modulation default) destroys symmetric 3-way
  crossings** — 100% below cutoff — because `p2` is *identically zero* for three
  equal orthogonal fibres. Realistic 3-fibre voxels survive but are dimmed 4x.
- **The kernel `p4` is structurally blind to crossings.**
  `Get_uniformly_distributed_SM_prior` draws `p4 = rand * p2 * 0.9`, so with
  `p2 = 0` the regression cannot represent `p4 = 0.54`.
- **Raw anisotropic power fails; noise does not cancel.** At SNR 15 it *inverts*:
  CSF 0.207 > orthogonal crossing 0.170.
- **The signal-domain route does not escape the noise floor.** Orthogonal
  crossing 0.0447 vs CSF 0.0365. Changing domain does not create information.
- **Lmax 8 is worse than Lmax 6** for the modulation work (unweighted FP 52.5%
  vs 21.9%). Note this is *not* true for deconvolution accuracy — section 6.4's
  Lmax 8 controls improve every arm. Different question, different answer.
- **`degenerate = 'clip'` is a bad modulation default.** In a blown-up voxel the
  raw `p` exceeds the clip so the voxel gets weight exactly 1.0, the maximum in
  the volume. `'reject'` exists. Still not changed.
- **`lambda_tikhonov` is inert.** Section 2.3, item 1. Stop sweeping it.

The one thing that has ever satisfied every constraint at once, still unshipped:
**noise-floor-subtracted anisotropic power over `l = 2,4`**, using `sigma` and
the kernel gains already in `out`. 0.0% FP on GM+CSF at both SNRs with all six
WM classes retained. Blocked on a free safety factor nobody has justified, and
it has never touched real data.

---

## 8. Next steps, roughly in order of value per hour

1. **Run `smi_manuscript_60deg.m` at full size.** Everything about it is built
   and checked, and it now runs all three arms, but **nothing has been run at
   `NREP = 1000`** — only at `SMOKE_TEST = true`, whose numbers are indicative
   only. The user runs it on their own hardware, which is faster than the
   container; do not burn hours running it here unless asked.

   What to read first: **Figure 6** (bias, spread and spurious peaks against
   SNR, one row per kernel and arm) and **Figure 7** (the three arms side by
   side as glyphs). The smoke run already reproduces the published ordering at
   60 degrees, Lmax 6, healthy kernel:

   | kernel | arm | SNR 10 correct / spurious | SNR 30 correct / spurious | SNR inf error |
   |---|---|---|---|---|
   | healthy | SMI | **100.0% / 0.000** | 100.0% / 0.000 | 1.72 deg |
   | healthy | SSST-CSD | 44.4% / 0.667 | 100.0% / 0.000 | **1.27 deg** (ceiling) |
   | healthy | ~~MSMT-CSD~~ | ~~85.2% / 0.000~~ | ~~92.6% / 0.000~~ | ~~7.80 deg~~ |
   | edema | SMI | 0.0% / 1.741 | **70.4% / 0.296** | 1.27 deg |
   | edema | SSST-CSD | 0.0% / 3.111 | 0.0% / 2.556 | 1.27 deg |
   | edema | ~~MSMT-CSD~~ | ~~33.3% / 0.074~~ | ~~48.1% / 0.000~~ | ~~11.72 deg~~ |

   **Every MSMT row above is void** -- it was measured with `-neg_lambda` at
   MRtrix's default of 1e-10, i.e. effectively unconstrained. See section 1.
   The SMI and SSST-CSD rows are unaffected. Re-run before quoting MSMT.

   Two things stand out, and the second is new. **In healthy tissue the arms
   trade places with SNR** — section 6.4's finding, now reproduced on one noise
   draw shared by all three arms instead of across separate runs. **In edema the
   trade disappears and SSST-CSD collapses**: 0% correct at both finite SNRs
   with 2.5-3 spurious peaks per voxel, where SMI still recovers 70% at SNR 30.
   That is the "less anisotropic signal, amplified in inverse proportion"
   mechanism (section 6.5) appearing as a measured method difference rather than
   as an argument. At `NREP = 27` these are 27 realisations — a smoke signal,
   not a result.

   **Its figures have still never been verified to RENDER at full size** — that
   gap already shipped one real bug (section 3, the shared-scale trap). The
   excuse is gone though: Octave's `print` works in this container (section 4).
   Render them before trusting them.
2. ~~**Wire up the CSD and MSMT-CSD arms.**~~ **Done, in the `.m` file itself.**
   Step 6b of `smi_manuscript_60deg.m` calls `dwi2fod csd` and
   `dwi2fod msmt_csd` on the same `S_noisy` the SMI arm is fitted to. The
   commented hook is gone and so is the two-language arrangement that briefly
   replaced it. See section 1.
3. **Revisit `examples/example_fODF_regularization_sweep.m`.** The user asked
   for this explicitly and it was deferred to its own patch: audit it for hidden
   bugs, and add 3D isometric panels of the reconstructed fODFs as
   `lambda_nonneg`, `lambda_tikhonov` and `tau` vary, in the style of the 2007
   CSD paper. Two things to settle with the user first: whether those panels
   sweep at one SNR or several, and whether they use the manuscript's 60 degree
   crossing or that script's existing 40/60/90 set.
4. **Regenerate the `Reports/` tables on the real HCP protocol.** Every number
   in `Reports/deconv_tables.md` and
   `Reports/REPORT_SMI_deconvolution_MonteCarlo.md` was measured on the older
   synthetic scheme, still on disk at `protocol/hcp_like_3shell.txt` and marked
   superseded. The pipeline now reads the real one, so this is a re-run rather
   than a code change.
2. **Run the drivers on real data and read the peak histogram.** Everything in
   every report is simulation. The single most informative number is
   `out.fODF_outlier.Ncap`: simulation says it should be **0** on a regularized
   fit, so a non-zero count means the regularization is not behaving on real
   data as it does in simulation, which is a bigger finding than the cap. If it
   is non-zero, look at `SMI_fODFcap_flagged.nii` — ventricles means the
   deconvolution is blowing up; brain edge means the pre-fit erosion is too
   gentle.
3. **Run `dwi2response dhollander` on a real brain** and compare against the
   phantom-derived responses in section 6.3. That is the one step of the Monte
   Carlo pipeline where "estimated from real data" is currently a synthetic
   stand-in, and it is the user's own stated methodology step 1.
4. **Test the MSMT shell-weighting hypothesis** (section 6.4). Refit MSMT-CSD
   with the low-b shells downweighted or dropped and see whether 45 degrees
   comes back. If it does, that is a genuine, publishable characterisation of
   when MSMT-CSD loses crossings, and it is maybe two hours of `dwi2fod` runs.
5. **Extend the Monte Carlo to unequal volume fractions and three-way
   crossings.** Everything measured so far is two equal fibres. Unequal
   fractions are where spurious-peak counts usually separate methods, and the
   three-way case is where the modulation work already knows SMI has trouble.
7. **Land `claude/freewater-simulations`**, or the measured numbers quoted in
   `Reports/REPORT_fODF_outlier_cap.md` have no reproducible source.
8. Spatial context beyond the cap — a neighbour-agreement *weight* rather than
   just an outlier test — is the remaining unexplored axis.

---

## 9. Documentation and clutter debt

- **`Reports/REPORT_fODF_modulation.md` §3 states a wrong cause** for the `pl4` loss.
  Two handoffs have now flagged it and nobody has fixed it. Fix it in place.
- **`Reports/REPORT_fODF_regularization_sweep.md` now disagrees with the shipped
  default.** It argues for 10; the default is 1. `SMI.m:935` and `README.md`
  both explain the disagreement, but the report itself does not, and a reader
  who finds it first will be misled. Add a header pointing at section 6 of
  `Reports/REPORT_SMI_deconvolution_MonteCarlo.md`.
- **Six reports, no index.** A reader has no map of which report supersedes
  which, or which numbers are simulation-only. Everything in every report is
  simulation-only; that is worth saying once, loudly, in `README.md`.
- **The patch files roughly double the apparent size of every change** because
  they duplicate their own diff. That is the repo convention, so it stays.
- **CI runs `claude-review` only.** There is no MATLAB or Octave runner, so
  `tests/test_fODF_outlier_cap.m`, `tests/test_SMI_outlier_cap.m` and
  `tests/test_SMI_response_helpers.m` are never executed by CI — they pass locally
  under Octave and that is invisible on a PR. An Octave job would be
  straightforward and the repo has never had one. With three test files now,
  this is worth more than it was.
- **Stale branches on origin:** `add-claude-github-actions-1785176603102`,
  `claude/github-repo-markdown-changes-t49gi0`,
  `claude/fodf-modulation-edema-38eh81` are all merged and can be deleted.

---

## 10. Working preferences observed

- Wants to be told when a suggestion is wrong, **with evidence**. Has pushed
  back correctly several times.
- **Prefers measured over assumed.** Every default in this repo that was
  guessed has turned out wrong by roughly 2x.
- **Distrusts results that look wrong, and says so.** When they did —
  *"I don't know if I fully trust your simulation setup, especially with the
  huge error in 45 degree crossings"* — they were half right, and the half they
  were right about was a real bug in a shipped default. Take that kind of
  pushback as a request to re-measure with a more trustworthy instrument, not as
  a request to re-explain.
- **Prefers established tools over reimplementations.** They explicitly asked
  for CSD and MSMT-CSD to go through MRtrix3 rather than a Python
  reimplementation, and that was the right call twice over (section 7).
  A mixed-language package is fine with them; the split they suggested is the
  one now in place.
- **Wants outsiders to be able to verify the work.** This is the driver behind
  the paused task in section 8, item 1, and it should shape how anything new is
  written.
- **Nothing goes in the repo until it is proven.** Stated explicitly.
- Asks where changes were made — point at file and line numbers.
- Wants patch files and README updates alongside every change.
- Figures matter; a manuscript is the eventual target. The manuscript-grade
  simulation is still blocked on ground-truth *density*: every voxel is a single
  kernel, so density-weighted scoring is meaningless until there is a second
  tissue type per voxel.
- Two driver scripts live only on the user's machine and arrive as file uploads.
  Peak truncation must run **before** modulation, and the outlier cap now
  supersedes it.

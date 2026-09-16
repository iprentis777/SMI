# Patch history

The raw development history of this fork, as mail-formatted patches, plus one
subfolder per feature that was built and later removed.

**New contributors do not need to read or apply any of this.** Start with the
root [`README.md`](../../README.md) and [`Additions/`](../../Additions/). This
directory exists so that every number in the reports has a provenance, and so
that a removed feature can be recovered rather than only described.

These patches duplicate changes already in Git history. Prefer `git log`.

## Feature removals

| folder | what | why it was removed |
|---|---|---|
| [`fODF_tikhonov/`](fODF_tikhonov/) | Tikhonov damping of the deconvolution | inert at the weights it shipped with, harmful above them |
| [`fODF_modulation/`](fODF_modulation/) | anisotropy modulation of the fODF | simulation-only evidence, unresolved convention questions |

Each holds the removal as reverse-applyable diffs, the measurements behind it,
and what would have to be proven before restoring it.

## Three things to know before applying anything

1. **Four patches are not `git am`-able.** `0019`–`0022` have no mail headers —
   they begin at `diff --git`. Use `git apply` for those four and
   `git am --3way` for the rest. An earlier version of this README said to use
   `git am --3way` throughout, which fails at `0019`.
2. **There are two `0010` patches**, and the numbers do not disambiguate them.
   Apply `0010simulationwalkthrough.patch` **before**
   `0010readmepipelineguidance.patch` — that is chronological order, though the
   timestamps look reversed because one is recorded in UTC and the other in
   PDT. The order column below is authoritative.
3. **`0019`–`0022` carry no commit message.** Every other patch explains the
   measurement behind its change; those four explain nothing, and `0022` is a
   52-file archival move. The summaries below are reconstructed from their
   diffs, not from the authors' own words.

## The patches, in application order

Categories: **method** changes the toolbox itself · **measurement** produces the
numbers a default rests on · **simulation** builds or refines the validation
harness · **figures** presentation correctness · **repo** organization and docs.

| # | patch | cat | what it does |
|---|---|---|---|
| 1 | `0001fODFregularization` | method | **Origin of both regularizers.** The deconvolution was unregularized least squares → large negative lobes at clinical SNR. Adds the non-negativity constraint *and* Tikhonov damping, both off by default. |
| 2 | `0002fODFregularizationsweep` | measurement | The sweep that turned "plausible round numbers" into measured ones, on synthetic two-fibre voxels with exactly known ground truth. |
| 3 | `0003fODFdefaultlambdanonneg` | method | Sets `lambda_nonneg = 10` from that sweep. **Later reverted to 1** — see the root README on why the two scores disagree. |
| 4 | `0004fODFmodulation` | method | Anisotropy modulation. A `p_00 = 1` fODF has an isotropic floor of 1/(4π) = 0.0796, above MRtrix's iFOD2 cutoff of 0.05, so 100% of simulated GM/CSF survives termination. **Since retired.** |
| 5 | `0005fODFoutliercap` | method | Post hoc cap on isolated pathologically bright glyphs. One-sided: never raises a voxel. |
| 6 | `0006deconvolutionMonteCarlo` | simulation | **Largest patch, 55 files.** The kernel-as-response bridge (`r_l = K_l(b)·√((2l+1)4π)` is literally an MRtrix response row) plus the original cross-language Monte Carlo pipeline. |
| 7 | `0007handoffrewrite` | repo | Agent handoff; records the CS_phase verification against MRtrix, correcting two earlier claims made against dipy. |
| 8 | `0008montecarloconfig` | simulation | Collapses two duplicate copies of "what the experiment is" into `mc_config.m`. Anti-drift. |
| 9 | `0009fileorganization` | repo | Creates `examples/`, `helpers/`, `tests/`, `Reports/`, `Figures/`. |
| 10 | `0010simulationwalkthrough` | simulation | The pipeline could be run but **not audited**. Adds a step-by-step notebook where every step ends in a `CHECK` against an independently computed value. |
| 11 | `0010readmepipelineguidance` | repo | README pipeline refocus. **Apply after the walkthrough** — see note 2. |
| 12 | `0011hcpprotocolandfigures` | simulation | Swaps the synthetic protocol for a **real HCP scheme** (288 volumes, 18 b≈0 + 90×3), exposing jittered b-values and b=5 s/mm² "b=0". |
| 13 | `0012manuscriptnotebook` | simulation | `smi_manuscript_60deg.m` — the walkthrough cut to manuscript configurations. |
| 14 | `0013manuscriptsnrsweep` | simulation | Single SNR=50 → swept `[5 10 20 30 50 Inf]`; NREP 25 → 1000. |
| 15 | `0014twokernelmanuscript` | simulation | Healthy and edema kernels in one run, sharing ground truth, protocol and noise seeds. |
| 16 | `0015figure1sharedscale` | figures | One shared radial scale so glyph *size* carries meaning. |
| 17 | `0016glyphaxislimits` | figures | 0015 had **no visible effect** — `axis equal` fixes aspect ratio, not limits, so autoscale silently undid it. |
| 18 | `0017scopescaleandhandoff` | figures | 0016 over-applied; scopes the pinned limits back to Figure 1 only. |
| 19 | `0018csdmsmtnotebook` | simulation | First working three-arm comparison (SMI / SSST-CSD / MSMT-CSD) driving **real MRtrix3 binaries**. |
| 20 | `0019mrtrixinmanuscript` | simulation | Moves the three arms out of the notebook into the `.m`. Adds `helpers/fODF_peak_score.m` so nothing can score an fODF a second way. *(no commit message)* |
| 21 | `0020msmtneglambdafix` | simulation | The `-neg_lambda` finding: `csd` and `msmt_csd` do not ship comparable defaults. Reclassifies "more noise must not help" as a note rather than a failure for bias-dominated arms. *(no commit message)* |
| 22 | `0021csdhealthyresponse` | simulation | CSD response built from the healthy kernel for both tissues; adds the response-derivation report. *(no commit message)* |
| 23 | `0022centraliseandglyphs` | repo | **52 operations, almost all renames** — moves the original cross-language pipeline to `Archive/deconv_pipeline/`. Archived, not deleted: those reports' numbers have no other provenance. *(no commit message)* |
| 24 | `0023glyphamplitudehandoff` | figures | Measures glyph spread over 32 noise realisations, and **reverts** code changes proposed during that investigation. |
| 25 | `0024humanreadablerepository` | repo | **227 insertions, 1019 deletions.** The previous outsider-friendliness pass; adds per-directory READMEs. |
| 26 | `0025whitemattersimulationrework` | simulation | Adds `smi_wm_60deg.m`, the fair two-arm design that is now the only active simulation. |

## Reading them as a story

Three sequences are worth knowing about:

- **`0001` → `0002` → `0003`:** a feature, the measurement that tunes it, and
  the default that measurement produced. This is the pattern the repository
  tries to hold to.
- **`0015` → `0016` → `0017`:** a fix, a fix to the fix, and a scoping of the
  second fix. `0015` changed a scale that `axis equal` silently discarded, and
  `0016`'s correction then over-applied. Worth reading together as a case of a
  change that looked right and did nothing.
- **`0004` and `0001`'s Tikhonov half** both ended up retired. Their measured
  outcomes live in the two subfolders above.

"""`dwi2response dhollander` -- unsupervised 3-tissue response estimation.

A reimplementation of MRtrix3's dhollander algorithm (Dhollander, Raffelt &
Connelly, ISMRM 2016, with the Dhollander, Mito, Raffelt & Connelly, ISMRM 2019
single-fibre WM refinement) following lib/mrtrix3/dwi2response/dhollander.py of
MRtrix 3.0.4 step by step. It is used here because the response function is the
one thing CSD and MSMT-CSD get that SMI does not, so how it is estimated has to
be the estimator people actually use, not a convenient stand-in.

DELIBERATE DEVIATIONS, all forced by running on a voxel population rather than
an image:

  * The brain mask is not eroded (`-erode 3`). The phantom has no brain edge,
    so there is nothing for the erosion to remove.
  * `mrthreshold`'s golden section search is replaced by an exact search over
    all candidate thresholds of the same cost function (mrtrix_ops.py).
  * The tensor is fitted by dipy rather than by `dwi2tensor`. Both are weighted
    linear least squares on the log signal over all shells.

The single-fibre WM step needs an MSMT-CSD fit with a deliberately sharp
"ewmrf" response, exactly as the MRtrix script does, and uses dipy's
MultiShellDeconvModel for it.

Also provides the two older selectors, so the comparison is not hostage to one
of them:

  fa        `dwi2response fa`, the top N voxels by FA
  tournier  `dwi2response tournier`, recursive calibration on the second-to-
            first peak ratio
"""
import sys
import numpy as np
from dipy.core.gradients import gradient_table
from dipy.core.sphere import HemiSphere
from dipy.data import get_sphere
from dipy.direction.peaks import peak_directions
from dipy.reconst import dti
from dipy.reconst.csdeconv import (AxSymShResponse,
                                   ConstrainedSphericalDeconvModel)
from dipy.reconst.mcsd import MultiShellDeconvModel, MultiShellResponse

import binio
from mrtrix_ops import optimal_threshold, threshold_top, amp2response

SH_CONST = 0.5 / np.sqrt(np.pi)          # Y_00, i.e. 1/sqrt(4 pi)
LMAX = 6
FA_THRESH = 0.2                          # dhollander -fa
PCT_SFWM = 0.5                           # dhollander -sfwm
PCT_GM = 2.0                             # dhollander -gm
PCT_CSF = 10.0                           # dhollander -csf


# ------------------------------------------------------------------ metrics
def signal_decay_metric(data, bvals):
    """The SDM of the dhollander script: a volume-weighted mean of
    log(mean b=0 / mean b) over the non-zero shells, capped at 10."""
    shells = np.unique(bvals)
    m0 = np.maximum(data[:, bvals == shells[0]], 0).mean(axis=1)
    num = np.zeros(data.shape[0])
    den = 0
    for b in shells[1:]:
        sel = bvals == b
        mb = np.maximum(data[:, sel], 0).mean(axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            sdm_b = np.log(m0 / mb)
        num += sdm_b * sel.sum()
        den += sel.sum()
    sdm = num / den
    ok = np.isfinite(sdm) & (sdm > 0)
    return np.where(ok, np.minimum(sdm, 10.0), 0.0), ok


def _thr_mask(values, submask):
    """`mrthreshold -mask submask`: threshold chosen from the masked voxels."""
    t = optimal_threshold(values[submask])
    out = np.zeros(values.shape, bool)
    out[submask] = values[submask] > t
    return out


# --------------------------------------------------------------- selection
def dhollander_voxels(data, bvals, bvecs, verbose=True):
    """Return dict of boolean masks: sfwm, gm, csf, plus the DTI eigenvectors.

    `data` is [Nvox x Ndwi], already divided by its own b=0 mean.
    """
    sdm, safe = signal_decay_metric(data, bvals)
    gtab = gradient_table(bvals * 1000.0, bvecs=bvecs)
    tf = dti.TensorModel(gtab).fit(data)
    fa = np.nan_to_num(tf.fa)
    vecs = tf.evecs[..., 0]

    # ---- crude segmentation
    crude_wm = safe & (fa > FA_THRESH)
    crude_nonwm = safe & ~crude_wm
    med_nonwm = np.median(sdm[crude_nonwm])
    crude_csf = _thr_mask(np.where(crude_nonwm, sdm - med_nonwm, 0.0), crude_nonwm)
    crude_gm = crude_nonwm & ~crude_csf

    # ---- refined segmentation
    med_wm = np.median(sdm[crude_wm])
    mad_wm = np.median(np.abs(sdm[crude_wm] - med_wm))
    outl = crude_wm & (sdm > med_wm + 1.4826 * mad_wm * 2.0)
    refined_wm = crude_wm & ~outl

    med_gm = np.median(sdm[crude_gm])
    gm_hi = crude_gm & (sdm > med_gm)
    gm_lo = crude_gm & ~gm_hi
    # in each half keep the voxels CLOSER to the median (mrthreshold -invert)
    hi_sel = np.zeros_like(gm_hi)
    if gm_hi.any():
        t = optimal_threshold((sdm - med_gm)[gm_hi])
        hi_sel[gm_hi] = (sdm - med_gm)[gm_hi] <= t
    lo_sel = np.zeros_like(gm_lo)
    if gm_lo.any():
        t = optimal_threshold((med_gm - sdm)[gm_lo])
        lo_sel[gm_lo] = (med_gm - sdm)[gm_lo] <= t
    refined_gm = np.where(gm_hi, hi_sel, lo_sel) & crude_gm

    csf_min = sdm[crude_csf].min() if crude_csf.any() else 0.0
    csf_extra = (outl & (sdm > csf_min)) | crude_csf
    refined_csf = _thr_mask(np.where(csf_extra, sdm - csf_min, 0.0), csf_extra)

    # ---- final voxel selection
    n_csf = int(round(refined_csf.sum() * PCT_CSF / 100.0))
    vox_csf = threshold_top(np.where(refined_csf, sdm, 0.0), n_csf)

    n_gm = int(round(refined_gm.sum() * PCT_GM / 100.0))
    med_rgm = np.median(sdm[refined_gm])
    vox_gm = threshold_top(
        np.where(refined_gm, np.abs(sdm - med_rgm) + 1.0, 0.0), n_gm, bottom=True)

    n_sfwm = int(round(refined_wm.sum() * PCT_SFWM / 100.0))
    vox_sfwm = _sfwm_2019(data, bvals, bvecs, refined_wm, vox_csf, n_sfwm)

    if verbose:
        print(f'  safe {safe.sum():6d}  crude WM {crude_wm.sum():6d}  '
              f'GM {crude_gm.sum():6d}  CSF {crude_csf.sum():6d}')
        print(f'  refined  WM {refined_wm.sum():6d}  GM {refined_gm.sum():6d}  '
              f'CSF {refined_csf.sum():6d}')
        print(f'  selected SFWM {vox_sfwm.sum():5d}  GM {vox_gm.sum():5d}  '
              f'CSF {vox_csf.sum():5d}')
    return {'sfwm': vox_sfwm, 'gm': vox_gm, 'csf': vox_csf, 'vecs': vecs,
            'fa': fa, 'sdm': sdm,
            'refined': {'wm': refined_wm, 'gm': refined_gm, 'csf': refined_csf}}


def _sfwm_2019(data, bvals, bvecs, refined_wm, vox_csf, n_sfwm):
    """The built-in Dhollander et al. (2019) single-fibre WM selection.

    Deconvolve the refined WM with a deliberately over-sharp WM response
    (`ewmrf`: alternating +/- coefficients, all of the same magnitude) plus the
    CSF response, and score each voxel by

        metric = peak amplitude of the WM FOD / (WM l=0 amplitude + CSF amplitude)

    which is high when the FOD is concentrated in one direction. Run at lmax 2
    to cut the candidates to 2N, then at lmax 6 to cut them to N.
    """
    from mtcsd import TwoTissueCSD
    from peaks import sh_basis

    shells = np.unique(bvals)
    mean_sig = data.mean(axis=1)
    coef = np.median(mean_sig[refined_wm]) * np.sqrt(4 * np.pi)
    r_csf = np.array([[np.maximum(data[vox_csf][:, bvals == b].mean(), 1e-9)
                       / SH_CONST] for b in shells])

    sphere = HemiSphere.from_sphere(get_sphere(name='repulsion724'))
    cand = np.where(refined_wm)[0]
    for lmax, keep in [(2, 2 * n_sfwm), (LMAX, n_sfwm)]:
        ncoef = len(range(0, lmax + 1, 2))
        # ewmrf.txt: the same magnitude at every order with alternating sign,
        # i.e. a deliberately over-sharp WM response
        ewm = np.zeros((len(shells), ncoef))
        for i, b in enumerate(shells):
            if b < 1e-6:
                ewm[i, 0] = coef
            else:
                ewm[i, :] = coef * np.array([(-1)**k for k in range(ncoef)])
        model = TwoTissueCSD(bvals, bvecs, ewm, [r_csf], lmax, sphere.vertices)
        Bs = sh_basis(sphere.vertices, lmax)
        nc = model.ncoef
        metric = np.zeros(len(cand))
        for j, v in enumerate(cand):
            x = model.fit(data[v])
            pk = float((Bs @ x[:nc]).max())
            denom = x[0] + x[nc]
            metric[j] = pk / denom if denom > 0 else 0.0
        order = np.argsort(metric)[::-1][:max(keep, 1)]
        cand = cand[order]
    out = np.zeros(refined_wm.shape, bool)
    out[cand] = True
    return out


# ----------------------------------------------------------- older selectors
def fa_voxels(data, bvals, bvecs, n_select=300):
    """`dwi2response fa`: the N voxels with the highest FA."""
    gtab = gradient_table(bvals * 1000.0, bvecs=bvecs)
    fa = np.nan_to_num(dti.TensorModel(gtab).fit(data).fa)
    m = threshold_top(fa, n_select)
    return m, dti.TensorModel(gtab).fit(data).evecs[..., 0]


def tournier_voxels(data, bvals, bvecs, n_select=300, n_iter=8, lmax=LMAX):
    """`dwi2response tournier`: recursive calibration on the peak ratio.

    Single shell (the highest b), as the command is. Each iteration fits CSD
    with the current response, scores every voxel by the ratio of its second
    peak to its first, keeps the N most single-fibre-like, and re-estimates.
    Returns (mask, peak directions of the kept voxels).
    """
    b = np.unique(bvals)
    keep = (bvals == 0) | (bvals == b[-1])
    gt = gradient_table(bvals[keep] * 1000.0, bvecs=bvecs[keep])
    sphere = HemiSphere.from_sphere(get_sphere(name='repulsion724'))
    d = data[:, keep]

    # MRtrix initialises from a prolate tensor of FA 0.7 -- here, from the
    # equivalent zonal response of a stick-like kernel
    from kernel import response_sh
    r_sh = response_sh(b[-1], [0.7, 2.0, 2.0, 0.2, 0.0], lmax)
    idx = np.arange(d.shape[0])
    dirs_keep = None
    for _ in range(n_iter):
        model = ConstrainedSphericalDeconvModel(
            gt, AxSymShResponse(1.0, r_sh), sh_order_max=lmax)
        ratio = np.ones(len(idx))
        first = np.zeros((len(idx), 3))
        for j, v in enumerate(idx):
            odf = model.fit(d[v]).odf(sphere)
            pd, pv, _ = peak_directions(odf, sphere,
                                        relative_peak_threshold=0.05,
                                        min_separation_angle=25)
            if len(pv) == 0:
                continue
            first[j] = pd[0]
            ratio[j] = pv[1] / pv[0] if len(pv) > 1 and pv[0] > 0 else 0.0
        order = np.argsort(ratio)[:n_select]
        new_idx, dirs_keep = idx[order], first[order]
        r_new = amp2response(d[new_idx][:, bvals[keep] == b[-1]],
                             bvecs[keep][bvals[keep] == b[-1]], dirs_keep, lmax)
        if np.allclose(r_new, r_sh, rtol=1e-3):
            r_sh, idx = r_new, new_idx
            break
        r_sh, idx = r_new, new_idx
    m = np.zeros(data.shape[0], bool)
    m[idx] = True
    return m, dirs_keep, idx


# ------------------------------------------------------------------- driver
def estimate(data, bvals, bvecs, algo='dhollander', verbose=True):
    """Return (response dict, selection dict). Responses are [Nshell x Ncoef]
    zonal harmonic arrays in MRtrix's convention, ready for MultiShellResponse.
    """
    shells = np.unique(bvals)
    ls = np.arange(0, LMAX + 1, 2)

    if algo == 'dhollander':
        sel = dhollander_voxels(data, bvals, bvecs, verbose=verbose)
        picks = {'wm': (sel['sfwm'], False), 'gm': (sel['gm'], True),
                 'csf': (sel['csf'], True)}
        vecs = sel['vecs']
    elif algo == 'fa':
        m, vecs = fa_voxels(data, bvals, bvecs)
        sel = {'sfwm': m, 'vecs': vecs}
        picks = {'wm': (m, False)}
    elif algo == 'tournier':
        m, dirs_keep, idx = tournier_voxels(data, bvals, bvecs)
        vecs = np.zeros((data.shape[0], 3))
        vecs[idx] = dirs_keep
        sel = {'sfwm': m, 'vecs': vecs}
        picks = {'wm': (m, False)}
    else:
        raise ValueError(algo)

    resp = {}
    for name, (mask, iso) in picks.items():
        R = np.zeros((len(shells), 1 if iso else len(ls)))
        for i, b in enumerate(shells):
            s = bvals == b
            if iso or b < 1e-6:
                v = amp2response(data[mask][:, s], bvecs[s], vecs[mask],
                                 LMAX, isotropic=True)
                R[i, 0] = v[0]
            else:
                R[i, :] = amp2response(data[mask][:, s], bvecs[s], vecs[mask], LMAX)
        resp[name] = R
    return resp, sel


def main(tag, algos=('dhollander', 'tournier', 'fa')):
    bvals = binio.load('bvals').ravel()
    bvecs = binio.load('bvecs')
    ph = binio.load('phantom_' + tag)
    label = binio.load('phantom_label_' + tag).ravel().astype(int)
    data = ph.reshape(-1, len(bvals), order='F')
    data = data / data[:, bvals == 0].mean(axis=1, keepdims=True)

    with open(binio.DATA + f'/phantom_classes_{tag}.txt') as f:
        names = [l.strip() for l in f if l.strip()]

    for algo in algos:
        print(f'\n=== dwi2response {algo} ' + '=' * 52)
        resp, sel = estimate(data, bvals, bvecs, algo=algo)
        m = sel['sfwm']
        print(f'  WM selection composition ({m.sum()} voxels):')
        for i, nm in enumerate(names, start=1):
            frac = 100 * np.sum(m & (label == i)) / max(m.sum(), 1)
            if frac > 0:
                print(f'     {nm:28s} {frac:6.1f}%   (volume '
                      f'{100*np.mean(label == i):5.1f}%)')
        for k, R in resp.items():
            binio.save(f'resp_{algo}_{k}_{tag}', R)
        r = resp['wm']
        print(f'  WM response, normalised to l=0 per shell:')
        for i in range(r.shape[0]):
            print('     ' + ' '.join(f'{v:9.4f}' for v in r[i] / r[i, 0]))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'p30',
         tuple(sys.argv[2:]) or ('dhollander', 'tournier', 'fa'))

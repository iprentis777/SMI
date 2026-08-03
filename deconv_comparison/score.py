"""Score every method on the Monte Carlo data.

Reads the SH coefficients each method wrote, maps SMI's into the same basis as
dipy's, runs the one peak finder in peaks.py over all of them, and reports

  * angular error of the peaks: median and mean (bias) and standard deviation
    (precision), both for the largest peak alone and averaged over all true
    fibres
  * the number of peaks found: missed fibres and spurious peaks
  * peak amplitude, normalised to each method's own noise-free median, so that
    methods whose fODFs are not on the same scale can still be compared
  * angular correlation coefficient against the band limited ground truth,
    which is scale free and uses the whole fODF rather than its maxima

BASIS. SMI's spherical harmonics and dipy's descoteaux07 basis span the same
space but are not the same basis. Rather than assume they agree, Octave writes
out its own SH matrix on the shared evaluation directions and the exact
[Nlm x Nlm] map between the two is solved for here. The residual of that solve
is printed: if it is not at round-off, everything downstream is suspect.
"""
import os
import sys
import numpy as np

import binio
import peaks as pk

LMAX = 6
METHODS = [('SMI constrained', 'sh_smiC_%s'),
           ('SMI unconstrained', 'sh_smiU_%s'),
           ('SSST-CSD', 'sh_csd_%s'),
           ('MSMT-CSD', 'sh_msmt_%s')]


def smi_to_dipy(tag):
    """[Nlm x Nlm] map from SMI's SH coefficients to dipy's, and its residual."""
    ev = binio.load('eval_dirs')
    Y_smi = binio.load('Y_smi_' + tag)
    B = pk.sh_basis(ev, LMAX)
    M, *_ = np.linalg.lstsq(B, Y_smi, rcond=None)
    resid = np.abs(B @ M - Y_smi).max()
    return M, resid


def load_all(tag, M):
    """Every method's fODF as dipy-basis SH coefficients, [Nvox x Nlm]."""
    out = {}
    for name, pat in METHODS:
        try:
            C = binio.load(pat % tag)
        except FileNotFoundError:
            continue
        if name.startswith('SMI'):
            C = C @ M.T
        out[name] = C
    return out


def summarise(tag, snr, base=None, verbose=True, cache=True):
    """Score one arm. Peak extraction over 40,000 voxels times four methods is
    a couple of minutes, and the tables and the figure both want it, so the
    per-voxel measurements are cached in data/ keyed by tag and method list.
    Delete data/score_<tag>.npz to force a recompute."""
    M, resid = smi_to_dipy(tag)
    if verbose:
        print(f'\nSMI -> dipy SH basis map: max|B M - Y| = {resid:.2e}')

    cache_file = os.path.join(binio.DATA, f'score_{tag}.npz')
    key = '|'.join(n for n, _ in METHODS)

    cond = binio.load('mc_cond_id_' + tag).ravel().astype(int)
    axes = binio.load('mc_gt_axes_' + tag)          # [2 x 3 x NCOND]
    angles = binio.load('mc_angles_' + tag).ravel()
    gt6 = binio.load('mc_sh_gt6_' + tag) @ M.T      # band limited truth
    ncond = len(angles)
    idx = {c: np.where(cond == c + 1)[0] for c in range(ncond)}
    true_ax = {c: [axes[k, :, c] for k in range(axes.shape[0])
                   if np.isfinite(axes[k, 0, c])] for c in range(ncond)}

    if cache and os.path.exists(cache_file):
        z = np.load(cache_file, allow_pickle=True)
        if str(z['key']) == key:
            res = z['res'].item()
            gt_err = z['gt_err'].item()
            gt_np = z['gt_np'].item()
            if verbose:
                _print(tag, snr, angles, idx, res, gt_err, gt_np, true_ax, base)
            return res, idx, angles, gt_err

    data = load_all(tag, M)
    B_sphere = pk.sh_basis(pk.SPHERE.vertices, LMAX)
    res = {}
    for name, C in data.items():
        P = pk.find_peaks(C, LMAX, B_sphere=B_sphere)
        npk = np.array([len(d) for d, _ in P])
        err = np.full(len(C), np.nan)     # mean error over ALL true fibres,
                                          # undefined when one was not found
        errp = np.full(len(C), np.nan)    # error of the LARGEST peak, always
                                          # defined as long as a peak exists
        amp = np.array([v[0] if len(v) else np.nan for _, v in P])
        for c in range(ncond):
            for v in idx[c]:
                d = P[v][0]
                _, e, ok = pk.match_error(d, true_ax[c])
                if ok:
                    err[v] = e
                if len(d):
                    errp[v] = min(pk.angle(d[0], a) for a in true_ax[c])
        a = pk.acc(C, gt6[cond - 1])
        res[name] = dict(npk=npk, err=err, errp=errp, amp=amp, acc=a)

    # the ground truth run through the identical pipeline: the ceiling
    Pg = pk.find_peaks(gt6, LMAX, B_sphere=B_sphere)
    gt_err, gt_np = {}, {}
    for c in range(ncond):
        d = Pg[c][0]
        gt_err[c] = (min(pk.angle(d[0], a) for a in true_ax[c])
                     if len(d) else np.nan)
        gt_np[c] = len(d)

    if cache:
        np.savez_compressed(cache_file, key=key, res=res, gt_err=gt_err,
                            gt_np=gt_np)
    if verbose:
        _print(tag, snr, angles, idx, res, gt_err, gt_np, true_ax, base)
    return res, idx, angles, gt_err


def _stats(x):
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan, np.nan
    return np.median(x), x.mean(), x.std()


def _print(tag, snr, angles, idx, res, gt_err, gt_np, true_ax, base):
    cn = ['single' if a == 0 else f'cross {a:g} deg' for a in angles]
    names = list(res.keys())
    w = 20

    print(f'\n{"="*(18+w*len(names))}')
    print(f'  SNR {snr}   {sum(len(v) for v in idx.values())} realisations')
    print(f'{"="*(18+w*len(names))}')

    print(f'\nangular error of the LARGEST peak, deg   (median [mean +/- sd]);'
          f'  "ceiling" is the band limited truth through the same peak finder')
    print(f'{"condition":>16s}' + ''.join(f'{n:>{w}s}' for n in names) +
          f'{"ceiling":>10s}')
    for c, nm in enumerate(cn):
        row = f'{nm:>16s}'
        for n in names:
            med, mu, sd = _stats(res[n]['errp'][idx[c]])
            row += f'{med:8.2f} [{mu:5.2f}+/-{sd:5.2f}]'
        row += f'{gt_err[c]:10.2f}'
        print(row)

    print(f'\nangular error averaged over ALL true fibres, deg  '
          f'(only the realisations where every fibre was found)')
    print(f'{"condition":>16s}' + ''.join(f'{n:>{w}s}' for n in names))
    for c, nm in enumerate(cn):
        row = f'{nm:>16s}'
        for n in names:
            med, mu, sd = _stats(res[n]['err'][idx[c]])
            row += f'{med:8.2f} [{mu:5.2f}+/-{sd:5.2f}]'
        print(row)

    print(f'\nfibres resolved: peaks found vs fibres present '
          f'(% of realisations)')
    print(f'{"condition":>16s}' + ''.join(f'{n:>{w}s}' for n in names) +
          f'{"ceiling":>10s}')
    for c, nm in enumerate(cn):
        ntrue = len(true_ax[c])
        row = f'{nm:>16s}'
        for n in names:
            k = res[n]['npk'][idx[c]]
            row += (f'{100*np.mean(k == ntrue):7.1f}% ok'
                    f'{100*np.mean(k > ntrue):6.1f}% +')
        row += f'{gt_np[c]:9d}'
        print(row)

    print(f'\nspurious peaks per realisation (mean number above the true count)')
    print(f'{"condition":>16s}' + ''.join(f'{n:>{w}s}' for n in names))
    for c, nm in enumerate(cn):
        ntrue = len(true_ax[c])
        row = f'{nm:>16s}'
        for n in names:
            k = res[n]['npk'][idx[c]]
            row += f'{np.mean(np.maximum(k - ntrue, 0)):{w}.3f}'
        print(row)

    print(f'\nangular correlation coefficient against the band limited truth '
          f'(mean +/- sd)')
    print(f'{"condition":>16s}' + ''.join(f'{n:>{w}s}' for n in names))
    for c, nm in enumerate(cn):
        row = f'{nm:>16s}'
        for n in names:
            a = res[n]['acc'][idx[c]]
            a = a[np.isfinite(a)]
            row += f'{a.mean():13.4f}+/-{a.std():5.3f}'
        print(row)

    if base is not None:
        print(f'\npeak amplitude relative to the same method noise free '
              f'(median [coefficient of variation])')
        print(f'{"condition":>16s}' + ''.join(f'{n:>{w}s}' for n in names))
        for c, nm in enumerate(cn):
            row = f'{nm:>16s}'
            for n in names:
                a = res[n]['amp'][idx[c]]
                a = a[np.isfinite(a)]
                b = base[n][c]
                row += f'{np.median(a)/b:13.3f} [{a.std()/max(np.mean(a),1e-12):4.2f}]'
            print(row)


def baseline(tag):
    """Median peak amplitude per method per condition on the noise free run."""
    res, idx, angles, _ = summarise(tag, 'noise free', verbose=False)
    return {n: [np.nanmedian(res[n]['amp'][idx[c]]) for c in range(len(angles))]
            for n in res}


if __name__ == '__main__':
    tag = sys.argv[1]
    snr = sys.argv[2] if len(sys.argv) > 2 else tag
    base = baseline(sys.argv[3]) if len(sys.argv) > 3 else None
    summarise(tag, snr, base=base)

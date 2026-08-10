"""Score the three deconvolutions, with MRtrix3's peak finder for all of them.

Every arm -- SMI, SSST-CSD, MSMT-CSD and the band limited ground truth -- is
written as an MRtrix SH image and run through the same `sh2peaks`, so peak
extraction is not something this package implements. What is implemented here
is only the bookkeeping: matching peaks to true fibre axes, counting them, and
the angular correlation coefficient.

ONE POST HOC RULE, applied identically to every arm. `sh2peaks -num 4` returns
up to four peaks whether or not they are meaningful, so a peak is kept only if
its ANISOTROPIC amplitude (its amplitude minus the fODF's isotropic term,
c_00*Y_00) is positive and at least REL_THRESH of the largest anisotropic
amplitude in that voxel. Subtracting the isotropic term is what makes the
threshold mean the same thing for SMI, whose fODF has unit mass and a fixed
1/(4*pi) floor, and for CSD, whose l=0 term is its apparent fibre density.
"""
import os
import sys
import numpy as np

import binio
import mrtrix_io as mio

MD = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mrtrix')
REL_THRESH = 0.30
Y00 = 0.5 / np.sqrt(np.pi)
ARMS = [('SMI constrained', 'smifod'),
        ('SSST-CSD', 'csdfod'),
        ('MSMT-CSD', 'msmtfod')]


def keep_peaks(dirs, amp, iso):
    """Apply the relative anisotropic-amplitude rule. Returns a boolean mask."""
    a = amp - iso[:, None]
    a = np.where(np.isfinite(a), a, -np.inf)
    top = a.max(axis=1, keepdims=True)
    with np.errstate(invalid='ignore'):
        keep = (a > 0) & (a >= REL_THRESH * np.maximum(top, 1e-30))
    return keep


def angle(u, v):
    return np.degrees(np.arccos(np.clip(np.abs(u @ v), 0, 1)))


def load_arm(stem, tag):
    sh = mio.read(os.path.join(MD, f'{stem}_{tag}.mih'))
    sh = sh.reshape(-1, sh.shape[-1], order='F')
    d, amp = mio.read_peaks(os.path.join(MD, f'pk_{stem}_{tag}.mih'))
    iso = sh[:, 0] * Y00
    keep = keep_peaks(d, amp, iso)
    return sh, d, amp, keep


def measure(tag):
    cond = binio.load('mc_cond_id_' + tag).ravel().astype(int)
    axes = binio.load('mc_gt_axes_' + tag)
    angles = binio.load('mc_angles_' + tag).ravel()
    ncond = len(angles)
    idx = {c: np.where(cond == c + 1)[0] for c in range(ncond)}
    true_ax = {c: [axes[k, :, c] for k in range(axes.shape[0])
                   if np.isfinite(axes[k, 0, c])] for c in range(ncond)}

    res = {}
    for name, stem in ARMS + [('ground truth', 'gtfod')]:
        sh, d, amp, keep = load_arm(stem, tag)
        n = len(sh)
        npk = keep.sum(axis=1)
        errp = np.full(n, np.nan)     # largest kept peak vs nearest true axis
        errm = np.full(n, np.nan)     # mean over all true axes, when all found
        pamp = np.full(n, np.nan)
        for c in range(ncond):
            ax = true_ax[c]
            for v in idx[c]:
                k = np.where(keep[v])[0]
                if k.size == 0:
                    continue
                order = k[np.argsort(-amp[v, k])]
                pamp[v] = amp[v, order[0]]
                errp[v] = min(angle(d[v, order[0]], a) for a in ax)
                if order.size >= len(ax):
                    used, errs = set(), []
                    for a in ax:
                        best, bi = np.inf, None
                        for i in order:
                            if i in used:
                                continue
                            e = angle(d[v, i], a)
                            if e < best:
                                best, bi = e, i
                        used.add(bi)
                        errs.append(best)
                    errm[v] = float(np.mean(errs))
        res[name] = dict(npk=npk, errp=errp, errm=errm, amp=pamp, sh=sh)

    # angular correlation coefficient against the band limited truth, l >= 2
    gt = res['ground truth']['sh']
    for name, _ in ARMS:
        # an arm fitted at a higher Lmax than the truth is compared on the
        # orders they share; the extra bands have no counterpart to correlate
        nc = min(res[name]['sh'].shape[1], gt.shape[1])
        u = res[name]['sh'][:, 1:nc]
        v = gt[:, 1:nc]
        num = (u * v).sum(1)
        den = np.sqrt((u * u).sum(1) * (v * v).sum(1))
        res[name]['acc'] = np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)
    return res, idx, angles, true_ax


def report(tag, snr):
    res, idx, angles, true_ax = measure(tag)
    cn = ['single' if a == 0 else f'cross {a:g} deg' for a in angles]
    names = [n for n, _ in ARMS]
    w = 20
    gt = res['ground truth']

    print(f'\n{"="*(18+w*len(names)+12)}')
    print(f'  SNR {snr}   {sum(len(v) for v in idx.values())} realisations'
          f'   (peaks by MRtrix sh2peaks)')
    print(f'{"="*(18+w*len(names)+12)}')

    print(f'\nangular error of the largest peak, deg (median [mean +/- sd])')
    print(f'{"condition":>16s}' + ''.join(f'{n:>{w}s}' for n in names) +
          f'{"ceiling":>10s}')
    for c, nm in enumerate(cn):
        row = f'{nm:>16s}'
        for n in names:
            x = res[n]['errp'][idx[c]]
            x = x[np.isfinite(x)]
            if x.size:
                row += f'{np.median(x):8.2f} [{x.mean():5.2f}+/-{x.std():5.2f}]'
            else:
                row += f'{"-":>20s}'
        g = gt['errp'][idx[c]]
        g = g[np.isfinite(g)]
        row += f'{np.median(g) if g.size else np.nan:10.2f}'
        print(row)

    print(f'\nfibres resolved (correct peak count), % of realisations')
    print(f'{"condition":>16s}' + ''.join(f'{n:>{w}s}' for n in names) +
          f'{"ceiling":>10s}')
    for c, nm in enumerate(cn):
        nt = len(true_ax[c])
        row = f'{nm:>16s}'
        for n in names:
            row += f'{100*np.mean(res[n]["npk"][idx[c]] == nt):{w}.1f}'
        row += f'{100*np.mean(gt["npk"][idx[c]] == nt):10.1f}'
        print(row)

    print(f'\nspurious peaks per realisation (mean above the true count)')
    print(f'{"condition":>16s}' + ''.join(f'{n:>{w}s}' for n in names))
    for c, nm in enumerate(cn):
        nt = len(true_ax[c])
        row = f'{nm:>16s}'
        for n in names:
            row += f'{np.mean(np.maximum(res[n]["npk"][idx[c]] - nt, 0)):{w}.3f}'
        print(row)

    print(f'\nangular correlation coefficient vs the band limited truth (mean +/- sd)')
    print(f'{"condition":>16s}' + ''.join(f'{n:>{w}s}' for n in names))
    for c, nm in enumerate(cn):
        row = f'{nm:>16s}'
        for n in names:
            a = res[n]['acc'][idx[c]]
            a = a[np.isfinite(a)]
            row += f'{a.mean():13.4f}+/-{a.std():5.3f}'
        print(row)
    return res, idx, angles


if __name__ == '__main__':
    report(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else sys.argv[1])

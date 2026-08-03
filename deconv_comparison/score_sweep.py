"""Score the non-negativity sweep, with MRtrix3's peak finder.

Same `sh2peaks` output, same keep-rule and same ground truth as
score_mrtrix.py, so these numbers sit on the same scale as the three-arm
comparison. Run after `sweep_nonneg.m` and `run_mrtrix.sh sweep <tag> <n>`.
"""
import os
import sys
import numpy as np

import binio
import score_mrtrix as sm


def main(tag, snr):
    with open(binio.DATA + f'/sweep_names_{tag}.txt') as f:
        names = [l.strip() for l in f if l.strip()]
    cond = binio.load('mc_cond_id_' + tag).ravel().astype(int)
    axes = binio.load('mc_gt_axes_' + tag)
    angles = binio.load('mc_angles_' + tag).ravel()
    ncond = len(angles)
    idx = {c: np.where(cond == c + 1)[0] for c in range(ncond)}
    true_ax = {c: [axes[k, :, c] for k in range(axes.shape[0])
                   if np.isfinite(axes[k, 0, c])] for c in range(ncond)}
    cn = ['single' if a == 0 else f'cross {a:g}' for a in angles]

    gt_sh, *_ = sm.load_arm('gtfod', tag)
    rows = {}
    for i, nm in enumerate(names, start=1):
        sh, d, amp, keep = sm.load_arm(f'sweep{i}fod', tag)
        n = len(sh)
        npk = keep.sum(axis=1)
        errp = np.full(n, np.nan)
        for c in range(ncond):
            for v in idx[c]:
                k = np.where(keep[v])[0]
                if k.size == 0:
                    continue
                best = k[np.argmax(amp[v, k])]
                errp[v] = min(sm.angle(d[v, best], a) for a in true_ax[c])
        u, w = sh[:, 1:], gt_sh[:, 1:]
        den = np.sqrt((u * u).sum(1) * (w * w).sum(1))
        acc = np.where(den > 0, (u * w).sum(1) / np.maximum(den, 1e-30), np.nan)
        rows[nm] = dict(npk=npk, errp=errp, acc=acc)

    print(f'\n{"="*96}\n  non-negativity sweep, SNR {snr}, '
          f'{len(cond)//ncond} realisations per condition, peaks by MRtrix '
          f'sh2peaks\n{"="*96}')
    for title, fn in [
            ('resolved (correct number of peaks), % of realisations',
             lambda r, c: f'{100*np.mean(r["npk"][idx[c]] == len(true_ax[c])):12.1f}%'),
            ('angular error of the largest peak, deg (median)',
             lambda r, c: f'{np.nanmedian(r["errp"][idx[c]]):13.2f}'),
            ('angular correlation coefficient vs the band limited truth (mean)',
             lambda r, c: f'{np.nanmean(r["acc"][idx[c]]):13.4f}'),
            ('spurious peaks per realisation (mean above the true count)',
             lambda r, c: f'{np.mean(np.maximum(r["npk"][idx[c]] - len(true_ax[c]), 0)):13.3f}')]:
        print(f'\n{title}')
        print(f'{"setting":24s}' + ''.join(f'{n:>13s}' for n in cn))
        for nm, r in rows.items():
            print(f'{nm:24s}' + ''.join(fn(r, c) for c in range(ncond)))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

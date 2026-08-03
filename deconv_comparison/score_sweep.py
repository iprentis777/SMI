"""Score the non-negativity weight sweep written by sweep_nonneg.m.

Same peak finder, same sphere and same ground truth as score.py, so the numbers
here sit on the same scale as the four-method comparison.
"""
import sys
import numpy as np

import binio
import peaks as pk
import score

LMAX = 6


def main(tag, snr):
    with open(binio.DATA + f'/sweep_names_{tag}.txt') as f:
        names = [l.strip() for l in f if l.strip()]
    M, resid = score.smi_to_dipy(tag)
    cond = binio.load('mc_cond_id_' + tag).ravel().astype(int)
    axes = binio.load('mc_gt_axes_' + tag)
    angles = binio.load('mc_angles_' + tag).ravel()
    gt6 = binio.load('mc_sh_gt6_' + tag) @ M.T
    ncond = len(angles)
    idx = {c: np.where(cond == c + 1)[0] for c in range(ncond)}
    true_ax = {c: [axes[k, :, c] for k in range(axes.shape[0])
                   if np.isfinite(axes[k, 0, c])] for c in range(ncond)}
    B_sphere = pk.sh_basis(pk.SPHERE.vertices, LMAX)
    cn = ['single' if a == 0 else f'cross {a:g}' for a in angles]

    rows = {}
    for i, nm in enumerate(names, start=1):
        C = binio.load(f'sh_sweep{i}_{tag}') @ M.T
        P = pk.find_peaks(C, LMAX, B_sphere=B_sphere)
        npk = np.array([len(d) for d, _ in P])
        errp = np.full(len(C), np.nan)
        for c in range(ncond):
            for v in idx[c]:
                d = P[v][0]
                if len(d):
                    errp[v] = min(pk.angle(d[0], a) for a in true_ax[c])
        a = pk.acc(C, gt6[cond - 1])
        rows[nm] = dict(npk=npk, errp=errp, acc=a)

    print(f'\n{"="*96}\n  non-negativity sweep, SNR {snr}, '
          f'{len(cond)//ncond} realisations per condition\n{"="*96}')
    print(f'\nresolved (correct number of peaks), % of realisations')
    print(f'{"setting":24s}' + ''.join(f'{n:>13s}' for n in cn))
    for nm, r in rows.items():
        row = f'{nm:24s}'
        for c in range(ncond):
            ntrue = len(true_ax[c])
            row += f'{100*np.mean(r["npk"][idx[c]] == ntrue):12.1f}%'
        print(row)

    print(f'\nangular error of the largest peak, deg (median)')
    print(f'{"setting":24s}' + ''.join(f'{n:>13s}' for n in cn))
    for nm, r in rows.items():
        row = f'{nm:24s}'
        for c in range(ncond):
            row += f'{np.nanmedian(r["errp"][idx[c]]):13.2f}'
        print(row)

    print(f'\nangular correlation coefficient vs the band limited truth (mean)')
    print(f'{"setting":24s}' + ''.join(f'{n:>13s}' for n in cn))
    for nm, r in rows.items():
        row = f'{nm:24s}'
        for c in range(ncond):
            row += f'{np.nanmean(r["acc"][idx[c]]):13.4f}'
        print(row)

    print(f'\nspurious peaks per realisation (mean above the true count)')
    print(f'{"setting":24s}' + ''.join(f'{n:>13s}' for n in cn))
    for nm, r in rows.items():
        row = f'{nm:24s}'
        for c in range(ncond):
            ntrue = len(true_ax[c])
            row += f'{np.mean(np.maximum(r["npk"][idx[c]] - ntrue, 0)):13.3f}'
        print(row)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

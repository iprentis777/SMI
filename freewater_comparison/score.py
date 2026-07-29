"""Scoring: SMI (regularized, unmodulated) vs CSD vs MSMT-CSD.

All methods are scored by the same code, on the same evaluation sphere, from
the same noisy data.
"""
import sys
import numpy as np
import binio
import common_score as cs

CUT_MRTRIX = 0.05           # MRtrix iFOD2 default -cutoff


def med_iqr(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan
    return np.median(x), np.percentile(x, 75) - np.percentile(x, 25)


def run(tag, snr_label):
    names = cs.cond_names()
    cid = cs.cond_id()
    axes = cs.gt_axes()
    NCOND = len(names)

    meth = {}
    for k, f in [('SMI', f'smi_amp_{tag}'),
                 ('CSD b3', f'csd_b3_amp_{tag}'),
                 ('CSD b1', f'csd_b1_amp_{tag}'),
                 ('MSMT', f'msmt_amp_{tag}')]:
        meth[k] = binio.load(f)

    idx = {c: np.where(cid == c + 1)[0] for c in range(NCOND)}
    ax = {c: [axes[k, :, c] for k in range(3) if np.isfinite(axes[k, 0, c])]
          for c in range(NCOND)}

    # ---- per voxel: peak error, peak count, peak amplitude
    res = {}
    for mk, M in meth.items():
        err = np.full(M.shape[1], np.nan)
        npk = np.zeros(M.shape[1], dtype=int)
        amp = np.zeros(M.shape[1])
        for c in range(NCOND):
            for v in idx[c]:
                col = M[:, v]
                amp[v] = cs.peak_amp(col)
                if ax[c]:
                    n, e, ok = cs.matched_error(col, ax[c])
                    npk[v], err[v] = n, (e if ok else np.nan)
                else:
                    npk[v] = len(cs.peaks_of(col)[0])
        res[mk] = dict(err=err, npk=npk, amp=amp)

    print('=' * 100)
    print(f'  SNR {snr_label}   -- SMI regularized (no modulation) vs CSD vs MSMT-CSD, all at Lmax 6')
    print('=' * 100)

    # ---------------------------------------------- A. angular accuracy
    print('\nA. PEAK ORIENTATION  median angular error, deg (IQR).  '
          '"nres" = reps that failed to resolve all fibres')
    print(f'{"condition":15s}' + ''.join(f'{m:>21s}' for m in meth))
    for c in range(NCOND):
        if not ax[c]:
            continue
        row = f'{names[c]:15s}'
        for mk in meth:
            e = res[mk]['err'][idx[c]]
            n = res[mk]['npk'][idx[c]]
            m, q = med_iqr(e)
            nres = int(np.sum(n < len(ax[c])))
            row += f'{m:9.2f} ({q:4.2f}) {nres:3d}' + ' ' * 2
        print(row)

    # ---------------------------------------------- B. free-water amplitude cost
    print('\nB. FREE WATER COST TO fODF AMPLITUDE   median peak(fw=0.4) / median peak(fw=0)')
    print(f'{"geometry":15s}' + ''.join(f'{m:>12s}' for m in meth))
    for g in ['single', 'cross60', 'cross45']:
        c0 = names.index(f'{g}_fw00')
        c4 = names.index(f'{g}_fw40')
        row = f'{g:15s}'
        for mk in meth:
            a0 = np.median(res[mk]['amp'][idx[c0]])
            a4 = np.median(res[mk]['amp'][idx[c4]])
            row += f'{a4/a0:12.3f}'
        print(row)

    # ---------------------------------------------- C. absolute amplitudes
    print('\nC. MEDIAN PEAK AMPLITUDE (raw fODF units, isotropic floor included)')
    print(f'{"condition":15s}' + ''.join(f'{m:>12s}' for m in meth))
    for c in range(NCOND):
        row = f'{names[c]:15s}'
        for mk in meth:
            row += f'{np.median(res[mk]["amp"][idx[c]]):12.4f}'
        print(row)

    # ---------------------------------------------- D. WM/CSF contrast
    ccsf = names.index('csf')
    print('\nD. TRACTOGRAPHY CONTRAST  (scale free)')
    print(f'{"":15s}{"WM+40%FW / CSF":>18s}{"CSF > 0.05":>13s}'
          f'{"adaptive cut":>14s}{"CSF above it":>14s}')
    print(f'{"":15s}{"peak ratio":>18s}{"(MRtrix dflt)":>13s}'
          f'{"keeps 95% WM":>14s}{"":>14s}')
    for mk in meth:
        a = res[mk]['amp']
        wm_fw = np.median(a[idx[names.index('single_fw40')]])
        csf = np.median(a[idx[ccsf]])
        # adaptive cutoff: retain 95% of ALL wm voxels incl. both fw levels
        wm_all = np.concatenate([a[idx[names.index(n)]] for n in names
                                 if n != 'csf'])
        T = np.percentile(wm_all, 5)
        fp_fixed = 100.0 * np.mean(a[idx[ccsf]] > CUT_MRTRIX)
        fp_adapt = 100.0 * np.mean(a[idx[ccsf]] > T)
        print(f'{mk:15s}{wm_fw/csf:18.2f}{fp_fixed:12.1f}%'
              f'{T:14.4f}{fp_adapt:13.1f}%')

    # ---------------------------------------------- E. compartment estimates
    print('\nE. FREE WATER FRACTION RECOVERY  (true 0.40 in the fw40 rows, '
          '0.00 in fw00, 0.95 in csf)')
    kern = binio.load(f'smi_kernel_{tag}')
    vf = binio.load(f'msmt_vf_{tag}')
    print(f'{"condition":15s}{"SMI fw":>16s}{"MSMT csf vf":>16s}'
          f'{"MSMT gm vf":>14s}{"MSMT wm vf":>14s}')
    for c in range(NCOND):
        m1, q1 = med_iqr(kern[idx[c], 4])
        print(f'{names[c]:15s}{m1:9.3f} ({q1:4.3f})'
              f'{np.median(vf[idx[c],0]):16.3f}{np.median(vf[idx[c],1]):14.3f}'
              f'{np.median(vf[idx[c],2]):14.3f}')
    print()


if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])

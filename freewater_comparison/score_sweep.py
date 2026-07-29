"""Tikhonov sweep: angular resolution against stability.

lambda_tikhonov was tuned for stability (REPORT_fODF_regularization_sweep.md)
and has never been traded off against angular resolution. This scores both
sides of that trade on identical data.
"""
import sys
import numpy as np
import binio
import common_score as cs


def run(tag, snr):
    names = cs.cond_names()
    cid = cs.cond_id()
    axes = cs.gt_axes()
    lam = binio.load(f'tik_lambdas_{tag}').ravel()
    idx = {c: np.where(cid == c + 1)[0] for c in range(len(names))}
    ax = {c: [axes[k, :, c] for k in range(3) if np.isfinite(axes[k, 0, c])]
          for c in range(len(names))}
    csf = names.index('csf')

    print(f'\n{"="*104}\n  TIKHONOV SWEEP, SNR {snr}   '
          f'(lambda_nonneg = 10 throughout, Lmax 6)\n{"="*104}')
    print(f'{"lambda":>8s}{"45cross resolved":>18s}{"45 err":>9s}'
          f'{"60 err":>9s}{"single err":>12s}{"CSF peak":>10s}'
          f'{"WM/CSF":>9s}{"CSF>0.05":>10s}')
    for il, L in enumerate(lam, start=1):
        M = binio.load(f'smi_amp_lam{il}_{tag}')
        amp = np.array([cs.peak_amp(M[:, v]) for v in range(M.shape[1])])

        def errs(cname):
            c = names.index(cname)
            e, n = [], 0
            for v in idx[c]:
                npk, er, ok = cs.matched_error(M[:, v], ax[c])
                if ok:
                    e.append(er)
                    n += 1
            return (np.median(e) if e else np.nan), n

        e45, n45 = errs('cross45_fw00')
        e60, _ = errs('cross60_fw00')
        e1, _ = errs('single_fw00')
        csf_pk = np.median(amp[idx[csf]])
        wm_all = np.concatenate([amp[idx[c]] for c in range(len(names))
                                 if c != csf])
        fp = 100 * np.mean(amp[idx[csf]] > 0.05)
        print(f'{L:8.2f}{n45:12d}/{len(idx[0]):<5d}'
              f'{e45:9.2f}{e60:9.2f}{e1:12.2f}{csf_pk:10.4f}'
              f'{np.median(amp[idx[names.index("single_fw40")]])/csf_pk:9.2f}'
              f'{fp:9.1f}%')


if __name__ == '__main__':
    for t, s in [('clean', 'noise free'), ('snr30', '30'), ('snr15', '15')]:
        if binio.exists(f'tik_lambdas_{t}'):
            run(t, s)

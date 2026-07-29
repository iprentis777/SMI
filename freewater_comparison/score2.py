"""Compartment-split edema: f held at 0.60, extra-axonal water converted to
free water step by step. Scores SMI against CSD and MSMT-CSD with both exact
and realistically-estimated response functions."""
import sys
import numpy as np
import binio
import common_score as cs

SUF = '_c'
METH = [('SMI', 'smi_amp_%s'),
        ('CSD exact', 'csd_b3_amp_%s'),
        ('CSD tournier', 'csdT_b3_amp_%s'),
        ('MSMT exact', 'msmt_amp_%s'),
        ('MSMT tensor', 'msmtE_amp_%s'),
        ('MSMT SH', 'msmtS_amp_%s')]


def run(tag, snr):
    names = cs.cond_names(SUF)
    cid = cs.cond_id(SUF)
    axes = cs.gt_axes(SUF)
    kern = binio.load('kern_gt' + SUF)
    idx = {c: np.where(cid == c + 1)[0] for c in range(len(names))}
    ax = {c: [axes[k, :, c] for k in range(3) if np.isfinite(axes[k, 0, c])]
          for c in range(len(names))}

    amp, err = {}, {}
    for mk, f in METH:
        M = binio.load(f % tag)
        amp[mk] = np.array([cs.peak_amp(M[:, v]) for v in range(M.shape[1])])
        e = np.full(M.shape[1], np.nan)
        for c in range(len(names)):
            if not ax[c]:
                continue
            for v in idx[c]:
                _, er, ok = cs.matched_error(M[:, v], ax[c])
                if ok:
                    e[v] = er
        err[mk] = e

    print(f'\n{"="*112}\n  COMPARTMENT SPLIT, SNR {snr}.  f = 0.60 fixed; '
          f'extra-axonal water converted to free water\n{"="*112}')

    for geo in ['single', 'cross60']:
        rows = [n for n in names if n.startswith(geo)]
        print(f'\n{geo}:  peak fODF amplitude, normalised to the fw = 0 voxel '
              f'of the same method')
        print(f'{"f_extra / fw":>14s}' + ''.join(f'{m:>15s}' for m, _ in METH))
        base = {mk: np.median(amp[mk][idx[names.index(rows[0])]]) for mk, _ in METH}
        for n in rows:
            c = names.index(n)
            fw = kern[idx[c][0], 4]
            lab = f'{1-0.6-fw:.2f} / {fw:.2f}'
            print(f'{lab:>14s}' + ''.join(
                f'{np.median(amp[mk][idx[c]])/base[mk]:15.3f}' for mk, _ in METH))

    print(f'\npeak angular error, deg (median)')
    print(f'{"condition":>16s}' + ''.join(f'{m:>15s}' for m, _ in METH))
    for n in names:
        if n == 'csf':
            continue
        c = names.index(n)
        print(f'{n:>16s}' + ''.join(
            f'{np.nanmedian(err[mk][idx[c]]):15.2f}' for mk, _ in METH))

    ccsf = names.index('csf')
    print(f'\ntractography contrast')
    print(f'{"":16s}{"edema/CSF":>12s}{"CSF peak":>11s}{"CSF>0.05":>11s}'
          f'{"CSF above 95%-WM cut":>22s}')
    ed = names.index('single_fe20')      # the f=.6 fw=.2 case
    for mk, _ in METH:
        a = amp[mk]
        wm_all = np.concatenate([a[idx[c]] for c in range(len(names)) if c != ccsf])
        T = np.percentile(wm_all, 5)
        print(f'{mk:16s}{np.median(a[idx[ed]])/np.median(a[idx[ccsf]]):12.2f}'
              f'{np.median(a[idx[ccsf]]):11.4f}'
              f'{100*np.mean(a[idx[ccsf]] > 0.05):10.1f}%'
              f'{100*np.mean(a[idx[ccsf]] > T):21.1f}%')

    # ---- weights: does an f-based weight dim these voxels at all?
    k = binio.load(f'smi_kernel_{tag}')
    vfE = binio.load(f'msmtE_vf_{tag}')
    print(f'\nrecovered compartment fractions (median)')
    print(f'{"f_extra / fw":>14s}{"true f":>9s}{"SMI f":>9s}{"SMI fw":>9s}'
          f'{"MSMT wm vf":>13s}{"MSMT csf vf":>13s}')
    for n in [x for x in names if x.startswith('single')] + ['csf']:
        c = names.index(n)
        fw = kern[idx[c][0], 4]
        lab = 'CSF' if n == 'csf' else f'{1-kern[idx[c][0],0]-fw:.2f} / {fw:.2f}'
        print(f'{lab:>14s}{kern[idx[c][0],0]:9.2f}{np.median(k[idx[c],0]):9.3f}'
              f'{np.median(k[idx[c],4]):9.3f}{np.median(vfE[idx[c],2]):13.3f}'
              f'{np.median(vfE[idx[c],0]):13.3f}')

    print(f'\nweight comparison (cutoff retains 95% of all WM incl. every '
          f'free-water level)')
    print(f'{"weight":24s}{"CSF surviving":>15s}{"edema(fw=.2)/healthy":>22s}'
          f'{"fw=.4 / healthy":>18s}{"worst WM kept":>15s}')
    pl = binio.load(f'smi_pl_{tag}')
    cands = {'none (unmodulated)': np.ones_like(amp['SMI']),
             'SMI  1-fw': np.clip(1 - k[:, 4], 0, 1),
             'SMI  f': np.clip(k[:, 0], 0, 1),
             'SMI  p2product': np.clip(k[:, 5] * pl[:, 0], 0, 1),
             'MSMT wm vol.frac.': np.clip(vfE[:, 2], 0, 1)}
    wm_c = [c for c in range(len(names)) if c != ccsf]
    h = names.index('single_fe40'); e2 = names.index('single_fe20')
    e4 = names.index('single_fe00')
    for wn, w in cands.items():
        a = amp['SMI'] * np.where(np.isfinite(w), w, 0)
        wm_all = np.concatenate([a[idx[c]] for c in wm_c])
        T = np.percentile(wm_all, 5)
        permed = np.array([np.median(a[idx[c]]) for c in wm_c])
        worst = wm_c[int(np.argmin(permed))]
        base = np.median(a[idx[h]])
        print(f'{wn:24s}{100*np.mean(a[idx[ccsf]] > T):14.1f}%'
              f'{np.median(a[idx[e2]])/base:22.2f}'
              f'{np.median(a[idx[e4]])/base:18.2f}'
              f'{100*np.mean(a[idx[worst]] > T):14.1f}%')


if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])

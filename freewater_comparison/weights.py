"""Extension: which per-voxel weight best separates CSF from WM-in-free-water?

Scored the way REPORT_fODF_modulation.md section 5.2 does it: pick the cutoff
that retains 95% of ALL white matter (both free water levels, so a weight
cannot win by discarding the edematous voxels -- that would force the cutoff
down), then measure what fraction of CSF survives it.

The weight is applied as a density modulation to SMI's own fODF, i.e.
amplitude -> w * amplitude, which is SMI.modulate_fODF's 'density' mode.
"""
import sys
import numpy as np
import binio
import common_score as cs


def run(tag, snr):
    names = cs.cond_names()
    cid = cs.cond_id()
    M = binio.load(f'smi_amp_{tag}')
    kern = binio.load(f'smi_kernel_{tag}')     # [f Da Depar Deperp fw p2 p4]
    pl = binio.load(f'smi_pl_{tag}')           # [pl2 pl4 pl6]
    vf = binio.load(f'msmt_vf_{tag}')          # [csf gm wm]

    amp = np.array([cs.peak_amp(M[:, v]) for v in range(M.shape[1])])
    idx = {c: np.where(cid == c + 1)[0] for c in range(len(names))}
    wm_c = [c for c, n in enumerate(names) if n != 'csf']
    ed_c = [c for c, n in enumerate(names) if n.endswith('fw40')]
    hl_c = [c for c, n in enumerate(names) if n.endswith('fw00')]
    csf = idx[names.index('csf')]

    def clip01(w):
        w = np.where(np.isfinite(w), w, 0.0)
        return np.clip(w, 0.0, 1.0)

    cands = {
        'none (unmodulated)': np.ones_like(amp),
        'SMI  1-fw': clip01(1 - kern[:, 4]),
        'SMI  f': clip01(kern[:, 0]),
        'SMI  p2product': clip01(kern[:, 5] * pl[:, 0]),
        'MSMT wm vol.frac.': clip01(vf[:, 2]),
    }

    print(f'\nF. WEIGHT COMPARISON, SNR {snr}   '
          f'(cutoff retains 95% of all WM incl. the 40% free water voxels)')
    print(f'{"weight":22s}{"CSF surviving":>15s}{"edema kept":>13s}'
          f'{"healthy kept":>14s}{"edema/healthy":>15s}')
    for k, w in cands.items():
        a = amp * w
        wm_all = np.concatenate([a[idx[c]] for c in wm_c])
        T = np.percentile(wm_all, 5)
        fp = 100 * np.mean(a[csf] > T)
        ed = 100 * np.mean(np.concatenate([a[idx[c]] for c in ed_c]) > T)
        hl = 100 * np.mean(np.concatenate([a[idx[c]] for c in hl_c]) > T)
        med_ed = np.median(np.concatenate([a[idx[c]] for c in ed_c]))
        med_hl = np.median(np.concatenate([a[idx[c]] for c in hl_c]))
        print(f'{k:22s}{fp:14.1f}%{ed:12.1f}%{hl:13.1f}%{med_ed/med_hl:15.2f}')


if __name__ == '__main__':
    for tag, snr in [('snr30', 30), ('snr15', 15)]:
        run(tag, snr)

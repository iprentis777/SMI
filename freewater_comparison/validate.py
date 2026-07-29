"""Noise-free validation. If the response construction or the conventions are
wrong, this is where it shows: with no noise every method should land on the
ground truth almost exactly."""
import numpy as np
import binio
import common_score as cs

names = cs.cond_names()
cid = cs.cond_id()
axes = cs.gt_axes()
gt = binio.load('gt_amp')                   # [Ndir x NCOND], noise free truth

# one representative voxel per condition (all reps identical when noise free)
rep = {c: int(np.where(cid == c + 1)[0][0]) for c in range(len(names))}

meth = {}
for k, f in [('SMI', 'smi_amp_clean'),
             ('CSD b3', 'csd_b3_amp_clean'),
             ('CSD b1', 'csd_b1_amp_clean'),
             ('MSMT', 'msmt_amp_clean')]:
    meth[k] = binio.load(f)
meth['GT L8'] = gt
meth['GT L6'] = binio.load('gt_amp6')

print('=' * 96)
print('NOISE-FREE VALIDATION  (Lmax 6 fit; ground truth generated at Lmax 8)')
print('=' * 96)
hdr = f'{"condition":16s} {"method":8s} {"npk":>4s} {"ang err":>9s} {"peak amp":>10s} {"sph mean":>10s}'
print(hdr)
print('-' * 96)
for c, nm in enumerate(names):
    ax = [axes[k, :, c] for k in range(3) if np.isfinite(axes[k, 0, c])]
    for mk, M in meth.items():
        col = M[:, c] if mk.startswith('GT') else M[:, rep[c]]
        npk, err, ok = cs.matched_error(col, ax) if ax else (len(cs.peaks_of(col)[0]), np.nan, True)
        print(f'{nm:16s} {mk:8s} {npk:4d} '
              f'{"" if np.isnan(err) else f"{err:7.2f}":>9s} '
              f'{cs.peak_amp(col):10.4f} {col.mean():10.4f}')
    print('-' * 96)

vf = binio.load('msmt_vf_clean')
print('\nMSMT volume fractions (noise free), per condition   [order: as returned]')
print(f'{"condition":16s} {"vf0":>8s} {"vf1":>8s} {"vf2":>8s}   true WM/FW')
for c, nm in enumerate(names):
    k = binio.load('kern_gt')[rep[c]]
    print(f'{nm:16s} {vf[rep[c],0]:8.3f} {vf[rep[c],1]:8.3f} {vf[rep[c],2]:8.3f}   '
          f'{1-k[4]:.2f} / {k[4]:.2f}')

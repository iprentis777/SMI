"""Figure: compartment-split edema, response estimation, and the sweeps."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dipy.reconst.shm import real_sh_descoteaux
from dipy.core.geometry import cart2sphere

import binio
import common_score as cs

SURFACE, INK, INK2, GRIDC = '#fcfcfb', '#0b0b0b', '#52514e', '#e3e2df'
C1, C2, C3, C4 = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
GT = '#8a8985'
plt.rcParams.update({
    'font.size': 8.5, 'axes.edgecolor': GRIDC, 'axes.labelcolor': INK2,
    'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
    'axes.facecolor': SURFACE, 'figure.facecolor': SURFACE, 'axes.grid': True,
    'grid.color': GRIDC, 'grid.linewidth': 0.6, 'axes.axisbelow': True})

SUF = '_c'
names = cs.cond_names(SUF)
cid = cs.cond_id(SUF)
axes_gt = cs.gt_axes(SUF)
idx = {c: np.where(cid == c + 1)[0] for c in range(len(names))}
FW = [0.0, 0.10, 0.20, 0.30, 0.40]

METH = [('SMI', 'smi_amp_%s', C1, '-'),
        ('CSD (b=3), estimated resp.', 'csdE_b3_amp_%s', C2, '-'),
        ('MSMT-CSD, exact resp.', 'msmt_amp_%s', C3, '-'),
        ('MSMT-CSD, estimated resp.', 'msmtE_amp_%s', C4, '-')]

fig, AX = plt.subplots(2, 2, figsize=(10.4, 7.6))
fig.subplots_adjust(hspace=0.46, wspace=0.29, top=0.815, bottom=0.09,
                    left=0.132, right=0.978)

# ---------------------------------------------- A: amplitude trajectory
ax = AX[0, 0]
for lab, f, col, ls in METH:
    M = binio.load(f % 'c_snr30')
    a = np.array([cs.peak_amp(M[:, v]) for v in range(M.shape[1])])
    y = [np.median(a[idx[names.index(f'single_fe{int(round((0.4-fw)*100)):02d}')]])
         for fw in FW]
    y = np.array(y) / y[0]
    ax.plot(FW, y, color=col, ls=ls, lw=2.0, marker='o', ms=5, label=lab)
ax.set_xlabel('free water fraction  (f held at 0.60; f_extra = 0.40 - fw)')
ax.set_ylabel('peak fODF amplitude, relative to fw = 0')
ax.set_ylim(0.7, 1.06)
ax.set_title('A   Edema as redistribution, not dilution', loc='left',
             fontsize=9.5, fontweight='bold', color=INK, pad=8)
ax.annotate('the case asked for', xy=(0.20, 0.905), xytext=(0.225, 0.775),
            fontsize=7, color=INK2,
            arrowprops=dict(arrowstyle='->', color=INK2, lw=0.8))

# ------------------------------------------------- B: angular error
ax = AX[0, 1]
for lab, f, col, ls in METH:
    M = binio.load(f % 'c_snr30')
    y = []
    for fw in FW:
        c = names.index(f'cross60_fe{int(round((0.4-fw)*100)):02d}')
        axl = [axes_gt[k, :, c] for k in range(3) if np.isfinite(axes_gt[k, 0, c])]
        e = [cs.matched_error(M[:, v], axl)[1] for v in idx[c]]
        y.append(np.nanmedian(np.array(e, float)))
    ax.plot(FW, y, color=col, ls=ls, lw=2.0, marker='o', ms=5, label=lab)
ax.set_xlabel('free water fraction')
ax.set_ylabel('60° crossing, median peak angular error (deg)')
ax.set_title('B   Response estimation costs MSMT 3-4x in angle',
             loc='left', fontsize=9.5, fontweight='bold', color=INK, pad=8)

# --------------------------- C: who lands in the response selection
ax = AX[1, 0]
comp = {'FA / healthy': [83.5, 16.5, 0.0, 0.0],
        'FA / edema': [86.5, 13.5, 0.0, 0.0],
        'anisotropy / healthy': [83.0, 16.5, 0.0, 0.5],
        'anisotropy / edema': [68.5, 14.0, 17.0, 0.5]}
labs = ['WM single', 'WM crossing', 'EDEMA', 'other']
cols = [C1, C3, C2, GRIDC]
y = np.arange(len(comp))
left = np.zeros(len(comp))
for k in range(4):
    v = np.array([comp[c][k] for c in comp])
    ax.barh(y, v, left=left, color=cols[k], height=0.62,
            label=labs[k], linewidth=0)
    for j, (l0, vv) in enumerate(zip(left, v)):
        if vv > 6:
            ax.text(l0 + vv/2, y[j], f'{vv:.0f}', ha='center', va='center',
                    fontsize=6.8, color='white' if k < 3 else INK2)
    left += v
ax.set_yticks(y)
ax.set_yticklabels(list(comp), fontsize=7.8)
ax.invert_yaxis()
ax.set_xlabel('composition of the top 5% "most anisotropic" voxels (%)')
ax.set_xlim(0, 100)
ax.grid(axis='y', b=False) if False else ax.yaxis.grid(False)
ax.legend(fontsize=6.8, frameon=False, ncol=4, loc='lower center',
          bbox_to_anchor=(0.5, -0.30), labelcolor=INK2, handlelength=1.2)
ax.set_title('C   Edema is invisible to an anisotropy selector',
             loc='left', fontsize=9.5, fontweight='bold', color=INK, pad=8)

# --------------------------------- D: the non-negativity sweep
ax = AX[1, 1]
nm2, cid2 = cs.cond_names(), cs.cond_id()
i2 = {c: np.where(cid2 == c + 1)[0] for c in range(len(nm2))}
c45 = nm2.index('cross45_fw00')
axl = [axes_gt2 for axes_gt2 in
       [cs.gt_axes()[k, :, c45] for k in range(3)]
       if np.isfinite(axes_gt2[0])]
lams = binio.load('nn_lambdas_clean').ravel()
res, errs = [], []
for il in range(1, len(lams) + 1):
    M = binio.load(f'smi_amp_nn{il}_clean')
    e = [cs.matched_error(M[:, v], axl) for v in i2[c45]]
    ok = [x[1] for x in e if x[2]]
    res.append(100 * len(ok) / len(i2[c45]))
    errs.append(np.median(ok) if ok else np.nan)
xs = np.arange(len(lams))
bars = ax.bar(xs, res, 0.55, color=[C1 if r > 50 else GRIDC for r in res],
              linewidth=0)
for x, r, e, L in zip(xs, res, errs, lams):
    lab = f'{e:.1f}°' if np.isfinite(e) else 'never\nresolved'
    ax.text(x, r + 3, lab, ha='center', fontsize=7,
            color=INK2 if np.isfinite(e) else '#e34948')
ax.set_xticks(xs)
ax.set_xticklabels([('off' if L == 0 else ('10\n(shipped)' if L == 10 else f'{L:g}')) for L in lams])
ax.set_xlabel('$\\lambda_{nonneg}$')
ax.set_ylabel('45° crossings resolved (%), NO NOISE')
ax.set_ylim(0, 118)
ax.set_title('D   Non-negativity closes the 45° crossing',
             loc='left', fontsize=9.5, fontweight='bold', color=INK, pad=8)

for a in (AX[0, 0], AX[0, 1]):
    a.legend(fontsize=6.8, frameon=False, labelcolor=INK2, handlelength=2.4)

fig.suptitle('Edema that spares the axons',
             x=0.078, y=0.968, ha='left', fontsize=13, fontweight='bold', color=INK)
fig.text(0.078, 0.878,
         'f = 0.60 held fixed while extra-axonal water is converted to free water, so no weight built on f can dim these voxels. '
         'SNR 30, 40 realisations,\nall methods at Lmax 6. "Estimated resp." = built from the top 5% most anisotropic voxels of a brain that contains 15% edema.',
         ha='left', fontsize=7.3, color=INK2)

fig.savefig('fodf_compartment_split.png', dpi=200, facecolor=SURFACE)
print('wrote fodf_compartment_split.png')

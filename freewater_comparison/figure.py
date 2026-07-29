"""Summary figure: SMI (regularized, unmodulated) vs CSD vs MSMT-CSD
recovering a WM fODF diluted with 40% free water."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

import binio
import common_score as cs

SURFACE = '#fcfcfb'
INK = '#0b0b0b'
INK2 = '#52514e'
GRID = '#e3e2df'
COL = {'SMI': '#2a78d6', 'CSD (b=3)': '#eb6834', 'MSMT-CSD': '#1baf7a'}
GT = '#8a8985'
FILES = {'SMI': 'smi_amp_%s', 'CSD (b=3)': 'csd_b3_amp_%s', 'MSMT-CSD': 'msmt_amp_%s'}

plt.rcParams.update({
    'font.size': 8.5, 'axes.edgecolor': GRID, 'axes.labelcolor': INK2,
    'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
    'axes.facecolor': SURFACE, 'figure.facecolor': SURFACE,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.axisbelow': True, 'xtick.direction': 'out', 'ytick.direction': 'out',
})

names = cs.cond_names()
cid = cs.cond_id()
axes_gt = cs.gt_axes()
idx = {c: np.where(cid == c + 1)[0] for c in range(len(names))}
ev = binio.load('eval_dirs')

amp = {}
for snr in ('snr30', 'snr15'):
    amp[snr] = {}
    for k, f in FILES.items():
        M = binio.load(f % snr)
        amp[snr][k] = np.array([cs.peak_amp(M[:, v]) for v in range(M.shape[1])])

fig, AX = plt.subplots(2, 2, figsize=(10.2, 7.4))
fig.subplots_adjust(hspace=0.44, wspace=0.28, top=0.835, bottom=0.09,
                    left=0.075, right=0.985)

# ------------------------------------------------------- A: fODF profile
ax = AX[0, 0]
n1 = axes_gt[0, :, 0]
th = np.degrees(np.arccos(np.clip(np.abs(ev @ n1), 0, 1)))
edges = np.arange(0, 91, 3.0)
ctr = 0.5 * (edges[:-1] + edges[1:])
bin_id = np.digitize(th, edges) - 1


def profile(M, cols):
    """Median-over-reps azimuthal average as a function of angle from n1."""
    p = np.full((len(ctr), len(cols)), np.nan)
    for j, v in enumerate(cols):
        a = M[:, v]
        for b in range(len(ctr)):
            m = bin_id == b
            if m.any():
                p[b, j] = a[m].mean()
    return np.nanmedian(p, axis=1)


c0, c4 = names.index('single_fw00'), names.index('single_fw40')
g6 = binio.load('gt_amp6')
gp = np.array([g6[bin_id == b, c0].mean() for b in range(len(ctr))])
ax.plot(ctr, gp, color=GT, lw=1.6, ls=(0, (1, 1.6)), zorder=1)
for k, f in FILES.items():
    M = binio.load(f % 'snr30')
    ax.plot(ctr, profile(M, idx[c0]), color=COL[k], lw=2.0)
    ax.plot(ctr, profile(M, idx[c4]), color=COL[k], lw=2.0, ls='--')
ax.set_xlabel('angle from true fibre axis (deg)')
ax.set_ylabel('fODF amplitude')
ax.set_xlim(0, 90)
ax.set_title('A   Single-fibre fODF profile, SNR 30', loc='left',
             fontsize=9.5, fontweight='bold', color=INK, pad=8)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], color=INK2, lw=1.8, label='no free water'),
                   Line2D([], [], color=INK2, lw=1.8, ls='--', label='40% free water'),
                   Line2D([], [], color=GT, lw=1.6, ls=(0, (1, 1.6)),
                          label='ground truth (Lmax 6)')],
          fontsize=7, frameon=False, loc='upper right', labelcolor=INK2,
          handlelength=2.4)

# --------------------------------------------- B: amplitude cost of water
ax = AX[0, 1]
geoms = ['single', 'cross60', 'cross45']
x = np.arange(len(geoms))
w = 0.24
for i, (k, f) in enumerate(FILES.items()):
    r = [np.median(amp['snr30'][k][idx[names.index(f'{g}_fw40')]]) /
         np.median(amp['snr30'][k][idx[names.index(f'{g}_fw00')]]) for g in geoms]
    pos = x + (i - 1) * (w + 0.02)
    ax.bar(pos, r, w, color=COL[k], label=k, linewidth=0)
    for p, v in zip(pos, r):
        ax.text(p, v + 0.02, f'{v:.2f}', ha='center', fontsize=6.6, color=INK2)
ax.axhline(1.0, color=INK2, lw=1.0, ls=':')
ax.axhline(0.60, color=INK2, lw=1.0, ls='--')
ax.text(-0.44, 1.02, 'no cost', fontsize=6.8, color=INK2, ha='left')
ax.text(-0.44, 0.505, 'cost = full tissue fraction (1 - fw = 0.60)',
        fontsize=6.8, color=INK2, ha='left')
ax.set_xticks(x)
ax.set_xticklabels(['single fibre', '60° crossing', '45° crossing'])
ax.set_ylabel('peak amplitude ratio,  40% FW / no FW')
ax.set_ylim(0, 1.18)
ax.set_title('B   What 40% free water costs the fODF amplitude', loc='left',
             fontsize=9.5, fontweight='bold', color=INK, pad=8)

# ------------------------------------------ C: tractography contrast
ax = AX[1, 0]
groups = [('healthy WM', 'single_fw00'), ('WM + 40% FW', 'single_fw40'),
          ('CSF', 'csf')]
x = np.arange(len(groups))
for i, (k, f) in enumerate(FILES.items()):
    v30 = [np.median(amp['snr30'][k][idx[names.index(c)]]) for _, c in groups]
    pos = x + (i - 1) * (w + 0.02)
    ax.bar(pos, v30, w, color=COL[k], label=k, linewidth=0)
    for p, v in zip(pos, v30):
        ax.text(p, v * 1.08, f'{v:.3f}', ha='center', fontsize=6.4, color=INK2)
ax.axhline(0.05, color='#e34948', lw=1.4, ls='--')
ax.text(-0.45, 0.061, 'MRtrix iFOD2 default cutoff 0.05', fontsize=6.8,
        color='#e34948', ha='left')
ax.set_yscale('log')
ax.set_ylim(0.012, 3.2)
ax.yaxis.set_major_locator(LogLocator(base=10, numticks=5))
ax.set_xticks(x)
ax.set_xticklabels([g for g, _ in groups])
ax.set_ylabel('median peak fODF amplitude (log)')
ax.set_title('C   Does the tracker stop where it should?  SNR 30', loc='left',
             fontsize=9.5, fontweight='bold', color=INK, pad=8)

# ---------------------------------------------- D: orientation accuracy
ax = AX[1, 1]
conds = [('single\nno FW', 'single_fw00'), ('single\n40% FW', 'single_fw40'),
         ('60° cross\nno FW', 'cross60_fw00'), ('60° cross\n40% FW', 'cross60_fw40')]
x = np.arange(len(conds))
for i, (k, f) in enumerate(FILES.items()):
    med, lo, hi = [], [], []
    M15 = binio.load(f % 'snr15')
    for _, c in conds:
        e = []
        for v in idx[names.index(c)]:
            ax_list = [axes_gt[j, :, names.index(c)] for j in range(3)
                       if np.isfinite(axes_gt[j, 0, names.index(c)])]
            _, err, ok = cs.matched_error(M15[:, v], ax_list)
            e.append(err if ok else np.nan)
        e = np.array(e, float)
        e = e[np.isfinite(e)]
        med.append(np.median(e)); lo.append(np.percentile(e, 25)); hi.append(np.percentile(e, 75))
    med, lo, hi = map(np.array, (med, lo, hi))
    pos = x + (i - 1) * (w + 0.02)
    ax.bar(pos, med, w, color=COL[k], label=k, linewidth=0,
           yerr=[med - lo, hi - med], error_kw=dict(ecolor=INK2, elinewidth=0.9,
                                                    capsize=2, capthick=0.9))
ax.set_xticks(x)
ax.set_xticklabels([g for g, _ in conds], fontsize=7.2)
ax.set_ylabel('peak angular error (deg), median & IQR')
ax.set_title('D   Orientation accuracy, SNR 15', loc='left',
             fontsize=9.5, fontweight='bold', color=INK, pad=8)

fig.suptitle('Recovering a white-matter fODF diluted with 40% free water',
             x=0.075, y=0.972, ha='left', fontsize=12.5, fontweight='bold', color=INK)
fig.text(0.075, 0.898,
         'HCP-like protocol (b = 0/1/2/3 ms·µm⁻², 18/90/90/90 dirs), Rician noise, 40 realisations per condition, all methods at Lmax 6.\n'
         'CSD and MSMT-CSD are given the exact ground-truth response functions; only SMI has to estimate its kernel from the data.',
         ha='left', fontsize=7.4, color=INK2)

from matplotlib.patches import Patch
fig.legend(handles=[Patch(facecolor=COL[k], label=k) for k in FILES],
           loc='upper right', bbox_to_anchor=(0.985, 0.978), ncol=3,
           frameon=False, fontsize=8.6, labelcolor=INK2, handlelength=1.5,
           columnspacing=1.4)

fig.savefig('fodf_freewater_comparison.png', dpi=200, facecolor=SURFACE)
print('wrote fodf_freewater_comparison.png')

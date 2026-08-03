"""Figure: the Monte Carlo comparison, as small multiples.

Rows are the three things measured -- angular accuracy, whether the right
number of fibres came back, and how often a peak appeared that is not there.
Columns are the fibre configuration. Each panel is scaled to its own condition,
because a single fibre and a 45 degree crossing differ by two orders of
magnitude in error and a shared row scale would flatten both; the y label names
the quantity and every panel in a row carries the same one. There is no second
y axis anywhere.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator

import score

COLORS = {'SMI constrained': '#2a78d6',
          'SMI unconstrained': '#eb6834',
          'SSST-CSD': '#1baf7a',
          'MSMT-CSD': '#eda100'}
INK = '#0b0b0b'
MUTED = '#52514e'
GRID = '#e6e5e2'


def style(ax):
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color('#c9c8c4')
    ax.tick_params(colors=MUTED, labelsize=8)


def collect(tags):
    """tags: list of (tag, snr_value). Returns nested dict of measurements."""
    out = {}
    angles = None
    ceiling = None
    for tag, snr in tags:
        res, idx, angles, gt_err = score.summarise(tag, snr, verbose=False)
        ceiling = gt_err
        for name, r in res.items():
            for c in range(len(angles)):
                sel = idx[c]
                e = r['errp'][sel]
                e = e[np.isfinite(e)]
                k = r['npk'][sel]
                a = r['acc'][sel]
                a = a[np.isfinite(a)]
                ntrue = 1 if angles[c] == 0 else 2
                out.setdefault(name, {}).setdefault(c, []).append(dict(
                    snr=snr,
                    err_med=np.median(e) if e.size else np.nan,
                    err_sd=e.std() if e.size else np.nan,
                    ok=100 * np.mean(k == ntrue),
                    spur=np.mean(np.maximum(k - ntrue, 0)),
                    acc=a.mean() if a.size else np.nan))
    return out, angles, ceiling


def main(out='fodf_deconv_montecarlo.png', *tagspec):
    tags = []
    for t in tagspec:
        tag, snr = t.split(':')
        tags.append((tag, float(snr)))
    data, angles, ceiling = collect(tags)
    names = list(data.keys())
    cn = ['single fibre' if a == 0 else f'crossing {a:g}°' for a in angles]
    ncol = len(angles)

    fig, axes = plt.subplots(3, ncol, figsize=(3.3 * ncol, 8.6), sharex=True)
    fig.patch.set_facecolor('#fcfcfb')
    rows = [('error of the largest peak (deg)', 'err_med', True),
            ('correct fibre count (%)', 'ok', False),
            ('spurious peaks per voxel', 'spur', False)]

    for ri, (ylab, key, logy) in enumerate(rows):
        for c in range(ncol):
            ax = axes[ri, c]
            for n in names:
                d = data[n][c]
                x = [r['snr'] for r in d]
                y = [r[key] for r in d]
                o = np.argsort(x)
                ax.plot(np.array(x)[o], np.array(y)[o], '-o', lw=2, ms=5,
                        color=COLORS.get(n, '#666'), label=n,
                        markeredgecolor='#fcfcfb', markeredgewidth=1.2)
            # the ceiling line is meaningless on a log axis when it is zero,
            # which it is for the single fibre: a Watson band limited to Lmax 6
            # still peaks exactly on its own axis, so there is nothing to draw
            if key == 'err_med' and np.isfinite(ceiling[c]) and ceiling[c] > 0.01:
                ax.axhline(ceiling[c], color=MUTED, lw=1.2, ls=':')
                ax.annotate('band-limited truth', (5, ceiling[c]),
                            fontsize=7, color=MUTED, xytext=(2, 3),
                            textcoords='offset points')
            if logy:
                ax.set_yscale('log')
            ax.set_xscale('log')
            ax.xaxis.set_major_locator(FixedLocator([5, 10, 20, 30, 50]))
            ax.xaxis.set_major_formatter(FixedFormatter(['5', '10', '20', '30', '50']))
            ax.xaxis.set_minor_locator(NullLocator())
            ax.set_xlim(4.4, 57)
            if ri == 0:
                ax.set_title(cn[c], fontsize=10, color=INK)
            if ri == 2:
                ax.set_xlabel('SNR', fontsize=9, color=MUTED)
            if c == 0:
                ax.set_ylabel(ylab, fontsize=9, color=MUTED)
            style(ax)

    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', ncol=len(names), frameon=False,
               fontsize=9, labelcolor=MUTED, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print('wrote', os.path.abspath(out))


if __name__ == '__main__':
    main(sys.argv[1], *sys.argv[2:])

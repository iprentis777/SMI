"""Figure: the SMI kernel as a response function, next to the estimated ones.

The glyphs and profiles here are the same quantities example_SMI_response_shview.m
draws in MATLAB. test_SMI_response_helpers.m checks that MATLAB side against
SMI's own forward model to 1e-15, and the response functions overlaid in panel F
are the ones `dwi2response` actually wrote.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import binio
from kernel import Kell, response_sh

LMAX = 8
KERN = [0.60, 2.0, 2.0, 0.50, 0.05]
SHELLS = [0.0, 1.0, 2.0, 3.0]
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
# l is an ordered variable, so it gets a sequential ramp rather than
# categorical hues; b is ordered too but only four shells are drawn and they
# are also identified by direct labels, so the categorical set is fine there.
SEQ = ['#0d3d70', '#1a5aa0', '#3a85d6', '#6ba3e0', '#9dc2ea']
INK = '#0b0b0b'
MUTED = '#52514e'


def zonal_profile(r, theta):
    x = np.cos(theta)
    P = {0: np.ones_like(x), 2: 1.5 * x**2 - .5,
         4: (35 * x**4 - 30 * x**2 + 3) / 8,
         6: (231 * x**6 - 315 * x**4 + 105 * x**2 - 5) / 16,
         8: (6435 * x**8 - 12012 * x**6 + 6930 * x**4 - 1260 * x**2 + 35) / 128}
    out = np.zeros_like(theta)
    for j, l in enumerate(range(0, 2 * (len(r) - 1) + 1, 2)):
        out += r[j] * np.sqrt((2 * l + 1) / (4 * np.pi)) * P[l]
    return out


def style(ax):
    ax.grid(True, color='#e6e5e2', lw=0.8)
    ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color('#c9c8c4')
    ax.tick_params(colors=MUTED, labelsize=8)


def main(tag='p30', out='fodf_response_shview.png'):
    th = np.linspace(0, np.pi, 361)
    R = np.stack([response_sh(b, KERN, LMAX) for b in SHELLS])

    fig = plt.figure(figsize=(12.5, 7.4))
    fig.patch.set_facecolor('#fcfcfb')

    # ---- A: rotational invariants
    ax = fig.add_subplot(2, 3, 1)
    bb = np.linspace(0, 3, 61)
    for j, l in enumerate(range(0, LMAX + 1, 2)):
        y = np.array([Kell(l, b, KERN)[0] for b in bb])
        ax.plot(bb, y, lw=2, color=SEQ[j])
        dy = {0: 0, 2: -3, 4: 7, 6: -4, 8: -12}[l]
        ax.annotate(f'$K_{l}$', (bb[-1], y[-1]), fontsize=9, color=MUTED,
                    xytext=(5, dy), textcoords='offset points')
    ax.set_xlim(0, 3.45)
    ax.axhline(0, color='#c9c8c4', lw=1)
    ax.set_xlabel('b  (ms/µm²)', fontsize=9, color=MUTED)
    ax.set_ylabel('$K_\\ell(b)$', fontsize=9, color=MUTED)
    ax.set_title('A   rotational invariants of the SMI kernel',
                 fontsize=10, color=INK, loc='left')
    style(ax)

    # ---- B: angular profile per shell
    ax = fig.add_subplot(2, 3, 2)
    for i, b in enumerate(SHELLS):
        y = zonal_profile(R[i], th)
        ax.plot(np.degrees(th), y, lw=2, color=C[i])
        ax.annotate(f'b = {b:g}', (np.degrees(th[40]), y[40]), fontsize=9,
                    color=MUTED, xytext=(4, 4), textcoords='offset points')
    ax.set_xlim(0, 180)
    ax.set_xlabel('angle from the fibre (deg)', fontsize=9, color=MUTED)
    ax.set_ylabel('R(θ)', fontsize=9, color=MUTED)
    ax.set_title('B   the same kernel as a response profile',
                 fontsize=10, color=INK, loc='left')
    style(ax)

    # ---- C: zonal coefficients, the MRtrix response rows
    ax = fig.add_subplot(2, 3, 3)
    ls = np.arange(0, LMAX + 1, 2)
    w = 0.2
    for i, b in enumerate(SHELLS):
        ax.bar(ls + (i - 1.5) * w, R[i], width=w, color=C[i],
               label=f'b = {b:g}', edgecolor='#fcfcfb', linewidth=0.8)
    ax.axhline(0, color='#c9c8c4', lw=1)
    ax.set_xticks(ls)
    ax.set_xticklabels([f'$r_{l}$' for l in ls])
    ax.set_ylabel('coefficient', fontsize=9, color=MUTED)
    ax.set_title('C   zonal harmonics = an MRtrix response file',
                 fontsize=10, color=INK, loc='left')
    ax.legend(frameon=False, fontsize=8, labelcolor=MUTED)
    style(ax)

    # ---- D: the shview glyphs
    axg = fig.add_subplot(2, 3, 4, projection='3d')
    tg = np.linspace(0, np.pi, 121)
    pg = np.linspace(0, 2 * np.pi, 121)
    TH, PH = np.meshgrid(tg, pg, indexing='ij')
    # each glyph is normalised to its own maximum, so the SHAPE is comparable
    # across shells; the amplitude drop with b is already panels B and C
    sep = 2.3
    for i in range(len(SHELLS)):
        A = zonal_profile(R[i], TH)
        A = A / np.abs(A).max()
        rr = np.abs(A)
        X = rr * np.sin(TH) * np.cos(PH) + (i - 1.5) * sep
        Y = rr * np.sin(TH) * np.sin(PH)
        Z = rr * np.cos(TH)
        norm = plt.Normalize(0, 1)
        axg.plot_surface(X, Y, Z, facecolors=plt.cm.Blues(0.30 + 0.6*norm(rr)),
                         rstride=2, cstride=2, linewidth=0, antialiased=True,
                         shade=True)
        axg.text((i - 1.5) * sep, 0, 1.25, f'b = {SHELLS[i]:g}',
                 ha='center', fontsize=8, color=MUTED)
    axg.set_box_aspect((4.6, 1.4, 1.4))
    axg.set_xlim(-3.6, 3.6); axg.set_ylim(-1.2, 1.2); axg.set_zlim(-1.2, 1.2)
    axg.view_init(elev=0, azim=-90)
    axg.set_axis_off()
    axg.set_title('D   response glyph per shell, as shview draws it',
                  fontsize=10, color=INK, loc='left')

    # ---- E: compartments
    ax = fig.add_subplot(2, 3, 5)
    f, Da, Dep, Dperp, fw = KERN
    parts = [('intra-axonal', [1, Da, 0, 0, 0], f),
             ('extra-axonal', [0, Da, Dep, Dperp, 0], 1 - f - fw),
             ('free water', [0, Da, Da, Da, 1], fw)]
    yt = zonal_profile(R[-1], th)
    ax.plot(np.degrees(th), yt, lw=2.4, color=INK)
    ax.annotate('total', (np.degrees(th[240]), yt[240]), fontsize=9, color=INK,
                xytext=(6, 8), textcoords='offset points')
    for j, (nm, kv, wgt) in enumerate(parts):
        y = wgt * zonal_profile(response_sh(SHELLS[-1], kv, LMAX), th)
        ax.plot(np.degrees(th), y, lw=1.8, ls='--', color=C[j])
        k = [120, 250, 300][j]
        ax.annotate(nm, (np.degrees(th[k]), y[k]), fontsize=8, color=MUTED,
                    xytext=[(-58, 6), (10, 10), (6, 12)][j],
                    textcoords='offset points')
    ax.set_ylim(-0.06, 0.85)
    ax.set_xlim(0, 180)
    ax.set_xlabel('angle from the fibre (deg)', fontsize=9, color=MUTED)
    ax.set_ylabel('R(θ)', fontsize=9, color=MUTED)
    ax.set_title('E   compartments of the response at b = 3',
                 fontsize=10, color=INK, loc='left')
    style(ax)

    # ---- F: SMI kernel vs the estimated responses
    ax = fig.add_subplot(2, 3, 6)
    # normalised by the l = 0 coefficient of each response, i.e. by its
    # spherical mean. Dividing by R(0) instead is meaningless here: at b = 3 the
    # response is nearly zero along the fibre.
    def norm_profile(r):
        return zonal_profile(r, th) * np.sqrt(4 * np.pi) / r[0]
    ax.plot(np.degrees(th), norm_profile(R[-1][:4]), lw=2.4, color=INK)
    ax.annotate('SMI kernel', (np.degrees(th[215]),
                norm_profile(R[-1][:4])[215]), fontsize=9, color=INK,
                xytext=(10, 4), textcoords='offset points')
    for j, algo in enumerate(['dhollander', 'tournier', 'fa']):
        try:
            r = binio.load(f'resp_{algo}_wm_{tag}')[-1].ravel()
        except FileNotFoundError:
            continue
        y = norm_profile(r)
        ax.plot(np.degrees(th), y, lw=1.8, color=C[j])
        k = [70, 95, 120][j]
        ax.annotate(f'{algo}', (np.degrees(th[k]), y[k]), fontsize=8,
                    color=MUTED, xytext=[(-56, -2), (-14, 10), (10, -2)][j],
                    textcoords='offset points')
    ax.set_xlim(0, 180)
    ax.set_xlabel('angle from the fibre (deg)', fontsize=9, color=MUTED)
    ax.set_ylabel('R(θ) normalised to its spherical mean', fontsize=9, color=MUTED)
    ax.set_title('F   estimated responses vs the SMI kernel, b = 3',
                 fontsize=10, color=INK, loc='left')
    style(ax)

    fig.tight_layout()
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print('wrote', os.path.abspath(out))


if __name__ == '__main__':
    import sys
    main(*sys.argv[1:])

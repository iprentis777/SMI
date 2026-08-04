"""Emit the result tables of REPORT_SMI_deconvolution_MonteCarlo.md as markdown.

The report quotes numbers; this writes them, so no number in it is transcribed
by hand. Peaks come from MRtrix's `sh2peaks` (score_mrtrix.py).

    python3 tables.py nf snr50:50 snr30:30 snr20:20 snr10:10 snr5:5
"""
import sys
import numpy as np

import score_mrtrix as sm


def fmt(x, n=2):
    return '—' if not np.isfinite(x) else f'{x:.{n}f}'


def main(nf_tag, *spec):
    tags = [(t.split(':')[0], float(t.split(':')[1])) for t in spec]
    snrs = ['noise free'] + [s for _, s in sorted(tags, key=lambda t: -t[1])]

    store = {}
    angles = ceiling = true_ax = None
    for tag, label in [(nf_tag, 'noise free')] + [(t, s) for t, s in tags]:
        res, idx, angles, true_ax = sm.measure(tag)
        store[label] = (res, idx)
    gtres, gtidx = store['noise free']
    base = {n: [np.nanmedian(gtres[n]['amp'][gtidx[c]])
                for c in range(len(angles))] for n, _ in sm.ARMS}
    ceiling = {c: np.nanmedian(gtres['ground truth']['errp'][gtidx[c]])
               for c in range(len(angles))}

    cn = ['single fibre' if a == 0 else f'crossing {a:g}°' for a in angles]
    names = [n for n, _ in sm.ARMS]

    def table(title, cell, ceil=False):
        print(f'\n**{title}**\n')
        for c, nm in enumerate(cn):
            extra = f' (band-limited truth: {fmt(ceiling[c])}°)' if ceil else ''
            print(f'\n_{nm}_{extra}')
            print('\n| SNR | ' + ' | '.join(names) + ' |')
            print('|---|' + '---|' * len(names))
            for s in snrs:
                res, idx = store[s]
                row = f'| {s} |'
                for n in names:
                    row += ' ' + cell(res[n], idx[c], c, n) + ' |'
                print(row)

    table('angular error of the largest peak, degrees (median [sd])',
          lambda r, sel, c, n: (f'{fmt(np.nanmedian(r["errp"][sel]))} '
                                f'[{fmt(np.nanstd(r["errp"][sel]))}]'), ceil=True)

    table('realisations returning the correct number of fibres, %',
          lambda r, sel, c, n: fmt(100 * np.mean(
              r['npk'][sel] == (1 if cn[c].startswith('single') else 2)), 1))

    table('spurious peaks per realisation (mean count above the truth)',
          lambda r, sel, c, n: fmt(np.mean(np.maximum(
              r['npk'][sel] - (1 if cn[c].startswith('single') else 2), 0)), 3))

    table('angular correlation coefficient vs the band limited truth (mean [sd])',
          lambda r, sel, c, n: (f'{fmt(np.nanmean(r["acc"][sel]), 4)} '
                                f'[{fmt(np.nanstd(r["acc"][sel]), 4)}]'))

    table('peak amplitude relative to the same method noise free '
          '(median [coefficient of variation])',
          lambda r, sel, c, n: (
              f'{fmt(np.nanmedian(r["amp"][sel]) / base[n][c], 3)} '
              f'[{fmt(np.nanstd(r["amp"][sel]) / max(np.nanmean(r["amp"][sel]), 1e-12), 2)}]'))


if __name__ == '__main__':
    main(sys.argv[1], *sys.argv[2:])

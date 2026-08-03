"""Emit the result tables of REPORT_SMI_deconvolution_MonteCarlo.md as markdown.

The report quotes numbers; this writes them, so no number in it is transcribed
by hand. Run after every SNR arm has been scored.

    python3 tables.py nf snr50:50 snr30:30 snr20:20 snr10:10 snr5:5
"""
import sys
import numpy as np

import score

ORDER = ['SMI constrained', 'SMI unconstrained', 'SSST-CSD', 'MSMT-CSD']


def fmt(x, n=2):
    return '—' if not np.isfinite(x) else f'{x:.{n}f}'


def main(nf_tag, *spec):
    tags = [(t.split(':')[0], float(t.split(':')[1])) for t in spec]
    base = score.baseline(nf_tag)

    all_res = {}
    angles = None
    ceiling = None
    for tag, snr in tags:
        res, idx, angles, gt = score.summarise(tag, snr, verbose=False)
        all_res[snr] = (res, idx)
        ceiling = gt
    nfres, nfidx, _, _ = score.summarise(nf_tag, 'nf', verbose=False)
    all_res['noise free'] = (nfres, nfidx)
    snrs = ['noise free'] + [s for _, s in sorted(tags, key=lambda t: -t[1])]

    cn = ['single fibre' if a == 0 else f'crossing {a:g}°' for a in angles]
    names = [n for n in ORDER if n in nfres]

    def table(title, cell, note=''):
        print(f'\n**{title}**{note}\n')
        for c, nm in enumerate(cn):
            print(f'\n_{nm}_ ' + (f'(band-limited truth: {fmt(ceiling[c])}°)'
                                  if 'error' in title else ''))
            print('\n| SNR | ' + ' | '.join(names) + ' |')
            print('|---|' + '---|' * len(names))
            for s in snrs:
                res, idx = all_res[s]
                row = f'| {s} |'
                for n in names:
                    row += ' ' + cell(res[n], idx[c], c, n) + ' |'
                print(row)

    table('angular error of the largest peak, degrees (median [sd])',
          lambda r, sel, c, n: (f'{fmt(np.nanmedian(r["errp"][sel]))} '
                                f'[{fmt(np.nanstd(r["errp"][sel]))}]'))

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

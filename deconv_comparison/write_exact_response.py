"""Write the EXACT response functions of the kernel that generated the data,
in MRtrix's response file format.

The control for section 7 of the report. A comparison in which MSMT-CSD comes
out worst on crossings is equally consistent with "the response it was given
was bad", so it is run a second time with the response no estimator could beat:
the delta-fODF response of the very kernel the signals were convolved with.

Rows are shells in ascending b, columns are the even zonal harmonic
coefficients r_l = K_l(b)*sqrt((2l+1)*4*pi), which is what `amp2response`
writes. The isotropic tissues get a single column, amplitude*sqrt(4*pi).
"""
import os
import sys
import numpy as np

from kernel import response_sh

LMAX = 6
K_WM = [0.60, 2.0, 2.0, 0.50, 0.02]     # the kernel gen_montecarlo.m uses
D_FW = 3.0
D_GM = 0.8
SHELLS = [0.0, 1.0, 2.0, 3.0]           # ms/um^2


def main(outdir):
    hdr = '# Shells: ' + ','.join(str(int(b * 1000)) for b in SHELLS) + '\n'
    wm = np.stack([response_sh(b, K_WM, LMAX) for b in SHELLS])
    gm = np.array([[np.exp(-b * D_GM) * np.sqrt(4 * np.pi)] for b in SHELLS])
    csf = np.array([[np.exp(-b * D_FW) * np.sqrt(4 * np.pi)] for b in SHELLS])
    for name, R in [('resp_wm_exact', wm), ('resp_gm_exact', gm),
                    ('resp_csf_exact', csf)]:
        with open(os.path.join(outdir, name + '.txt'), 'w') as f:
            f.write(hdr)
            for row in R:
                f.write(' '.join(f'{v:.12g}' for v in row) + '\n')
        print(f'{name}.txt')
    # the single-shell response for dwi2fod csd is the b = 3 row alone
    with open(os.path.join(outdir, 'resp_wm_exact_b3.txt'), 'w') as f:
        f.write('# Shells: 3000\n')
        f.write(' '.join(f'{v:.12g}' for v in wm[-1]) + '\n')
    print('resp_wm_exact_b3.txt')
    print('WM exact, normalised to l=0 per shell:')
    for i, b in enumerate(SHELLS):
        print(f'  b={b:g}  ' + ' '.join(f'{v:8.4f}' for v in wm[i] / wm[i, 0]))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'mrtrix')

"""Estimate CSD / MSMT-CSD response functions the way the real tools do.

Both take their WM response from the most anisotropic voxels in the volume.
That is normally safe: the most anisotropic voxels are coherent single-fibre
white matter. This script asks whether it is still safe when edema preserves
the intra-axonal fraction and only converts extra-axonal water into free water.

Selection is done explicitly rather than through dipy's threshold heuristics,
which are tuned for real brains and mis-fire on a synthetic volume (they leave
the GM mask empty here). The WM rule is the literal one -- top 5% by FA.
dipy's own `mask_for_response_msmt` is reported alongside for reference.
"""
import sys
import numpy as np
from dipy.core.gradients import gradient_table
from dipy.reconst.csdeconv import recursive_response
from dipy.reconst.mcsd import response_from_mask_msmt
import dipy.reconst.dti as dti

import binio
from kernel import response_sh

LMAX = 6
WM_PCT = 95.0          # "most anisotropic": top 5% by FA
CLASSES = {
    'healthy': ['WM single', 'WM crossing', 'WM partial vol', 'GM', 'CSF'],
    'edema': ['WM single', 'WM crossing', 'EDEMA f=.6 fw=.2',
              'WM partial vol', 'GM', 'CSF'],
}


def main(tag, variant):
    bvals = binio.load('bvals').ravel()
    bvecs = binio.load('bvecs')
    brain = binio.load(f'brain_{tag}')
    label = binio.load(f'brain_label_{tag}').ravel().astype(int)
    shape = brain.shape[:3]
    names = CLASSES[variant]

    b0 = brain[..., bvals == 0].mean(axis=-1, keepdims=True)
    data = brain / b0
    gtab = gradient_table(bvals * 1000.0, bvecs=bvecs)

    # DTI for the selection is fitted on b <= 1 only; a tensor fitted through
    # b = 3 has badly inflated FA and leaves no voxel below the GM threshold.
    low = (bvals == 0) | (bvals == 1.0)
    gt_low = gradient_table(bvals[low] * 1000.0, bvecs=bvecs[low])
    tf = dti.TensorModel(gt_low).fit(data[..., low])
    fa = tf.fa.reshape(-1, order='F')
    md = tf.md.reshape(-1, order='F')

    m_wm = fa >= np.percentile(fa, WM_PCT)
    m_gm = (fa < 0.15) & (md < 1.5e-3)
    m_csf = (fa < 0.20) & (md > 2.0e-3)

    print(f'\n{"="*88}\n  {variant.upper()} BRAIN -- what lands in each response mask\n{"="*88}')
    print(f'{"mask":26s}{"n":>6s}   ' + ''.join(f'{n:>19s}' for n in names))
    for k, msk in [('WM (top 5% FA)', m_wm), ('GM', m_gm), ('CSF', m_csf)]:
        n = int(msk.sum())
        row = f'{k:26s}{n:6d}   '
        for i in range(1, len(names) + 1):
            row += f'{100*np.sum(msk & (label == i))/max(n,1):18.1f}%'
        print(row)
    print(f'{"share of whole volume":26s}{len(label):6d}   ' +
          ''.join(f'{100*np.mean(label == i):18.1f}%'
                  for i in range(1, len(names) + 1)))

    def as3(m):
        return m.reshape(shape, order='F')

    # ------------------------------------------------- CSD, single shell b=3
    keep = (bvals == 0) | (bvals == 3.0)
    gt3 = gradient_table(bvals[keep] * 1000.0, bvecs=bvecs[keep])
    rr = recursive_response(gt3, data[..., keep], mask=as3(m_wm),
                            sh_order_max=LMAX, peak_thr=0.01, init_fa=0.08,
                            init_trace=0.0021, iter=8, convergence=0.001,
                            parallel=False)
    r_csd = np.asarray(rr.dwi_response, dtype=float)
    binio.save(f'resp_csd_{tag}', r_csd)

    # -------------------------------------------------------- MSMT, 3 tissue
    wm_rf, gm_rf, csf_rf = response_from_mask_msmt(
        gtab, data, as3(m_wm), as3(m_gm), as3(m_csf))
    for nm, arr in [('wm', wm_rf), ('gm', gm_rf), ('csf', csf_rf)]:
        binio.save(f'resp_msmt_{nm}_{tag}', np.asarray(arr, dtype=float))

    TISSUE = [0.60, 2.0, 2.0, 0.50, 0.0]
    r_true = response_sh(3.0, TISSUE, LMAX)
    print(f'\nCSD response at b=3, normalised to l=0:')
    print('  l              ' + ''.join(f'{l:>10d}' for l in range(0, LMAX + 1, 2)))
    print('  exact (fw=0)   ' + ''.join(f'{v:10.4f}' for v in r_true / r_true[0]))
    print('  estimated      ' + ''.join(f'{v:10.4f}' for v in r_csd / r_csd[0]))
    w = np.asarray(wm_rf)
    print(f'\nMSMT WM tensor evals x1e3, per shell (b = 1/2/3):')
    for i in range(w.shape[0]):
        print(f'  shell {i}: {np.round(w[i, :3]*1e3, 4)}   S0 {w[i, 3]:.4f}')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

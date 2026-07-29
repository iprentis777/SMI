"""Response estimation from "the most anisotropic voxels", two ways.

Selector 'fa'    : top 5% by tensor FA, fitted on b <= 1. This is what dipy's
                   `auto_response_ssst` and `mask_for_response_msmt` do.
Selector 'aniso' : top 5% by the anisotropy of the b = 3 shell signal itself,
                   ||S_2|| / S_0. This is the high-b analogue and is the one
                   the free water cannot hide from, since free water is gone by
                   b = 3 and only affects S_0 through the b = 0 normalisation.

Both then estimate the response the way every tool does: rotate each selected
voxel onto a common axis and average, per shell, per l.

A peak-ratio selector (`tournier.py`, the `dwi2response tournier` criterion) was
tried and is NOT used here: at Lmax 6 a 50 degree crossing does not resolve, so
it presents a single peak, scores a second-to-first ratio of exactly 0, and is
selected as the MOST single-fibre-like voxel in the volume. In the healthy brain
it filled 88% of the selection with crossings. That degeneracy is a real
property of the criterion at this angular order, but it makes the resulting
response useless as a comparison point.
"""
import sys
import numpy as np
from dipy.core.geometry import cart2sphere, vec2vec_rotmat
from dipy.core.gradients import gradient_table
from dipy.reconst.shm import (real_sh_descoteaux, real_sh_descoteaux_from_index)
import dipy.reconst.dti as dti

import binio
from kernel import response_sh

SH_CONST = 0.5 / np.sqrt(np.pi)
LMAX = 6
PCT = 95.0
CLASSES = {
    'healthy30': ['WM single', 'WM crossing', 'WM partial vol', 'GM', 'CSF'],
    'edema30': ['WM single', 'WM crossing', 'EDEMA f=.6 fw=.2',
                'WM partial vol', 'GM', 'CSF'],
}


def shell_anisotropy(data, bvals, bvecs, b):
    sel = bvals == b
    _, th, ph = cart2sphere(*bvecs[sel].T)
    B, _, _ = real_sh_descoteaux(4, th, ph)
    c = np.linalg.lstsq(B, data[:, sel].T, rcond=None)[0].T
    return np.linalg.norm(c[:, 1:6], axis=1) / np.abs(c[:, 0])


def main(tag):
    bvals = binio.load('bvals').ravel()
    bvecs = binio.load('bvecs')
    brain = binio.load(f'brain_{tag}')
    label = binio.load(f'brain_label_{tag}').ravel().astype(int)
    names = CLASSES[tag]
    b0 = brain[..., bvals == 0].mean(axis=-1, keepdims=True)
    data = (brain / b0).reshape(-1, len(bvals), order='F')

    low = (bvals == 0) | (bvals == 1.0)
    tf = dti.TensorModel(gradient_table(bvals[low] * 1000.0,
                                        bvecs=bvecs[low])).fit(data[:, low])
    fa, md = tf.fa, tf.md
    pev = tf.evecs[..., 0]
    a3 = shell_anisotropy(data, bvals, bvecs, 3.0)

    sels = {'fa': fa >= np.percentile(fa, PCT),
            'aniso': a3 >= np.percentile(a3, PCT)}

    print(f'\n{"="*94}\n  {tag}: who lands in the top 5% "most anisotropic"\n{"="*94}')
    print(f'{"selector":10s}{"n":>6s}   ' + ''.join(f'{n:>17s}' for n in names))
    print(f'{"(volume)":10s}{len(label):6d}   ' +
          ''.join(f'{100*np.mean(label==i):16.1f}%'
                  for i in range(1, len(names) + 1)))
    for k, m in sels.items():
        n = int(m.sum())
        print(f'{k:10s}{n:6d}   ' + ''.join(
            f'{100*np.sum(m & (label == i))/n:16.1f}%'
            for i in range(1, len(names) + 1)))

    m_gm = (fa < 0.15) & (md < 1.5e-3)
    m_csf = (fa < 0.20) & (md > 2.0e-3)
    shells = np.unique(bvals)
    ls = np.arange(0, LMAX + 1, 2)

    for k, m in sels.items():
        sel_idx = np.where(m)[0]
        R = np.zeros((len(shells), 2 + len(ls)))
        for i, b in enumerate(shells):
            s = bvals == b
            if b == 0:
                R[i, 2] = data[sel_idx][:, s].mean() / SH_CONST
            else:
                acc = np.zeros(len(ls))
                for j in sel_idx:
                    Rm = vec2vec_rotmat(pev[j], np.array([0, 0, 1.0]))
                    g = (Rm @ bvecs[s].T).T
                    _, th, ph = cart2sphere(*g.T)
                    B = real_sh_descoteaux_from_index(0, ls, th[:, None], ph[:, None])
                    acc += np.linalg.lstsq(B, data[j, s], rcond=-1)[0]
                R[i, 2:] = acc / len(sel_idx)
            R[i, 0] = data[m_csf][:, s].mean() / SH_CONST
            R[i, 1] = data[m_gm][:, s].mean() / SH_CONST
        binio.save(f'resp_{k}_{tag}', R)

    print(f'\nWM response at b=3, normalised to l=0:')
    ex = response_sh(3.0, [0.60, 2.0, 2.0, 0.50, 0.0], LMAX)
    print('  exact (delta fODF) ' + ''.join(f'{v:9.3f}' for v in ex / ex[0]))
    for k in sels:
        R = binio.load(f'resp_{k}_{tag}')
        print(f'  estimated ({k:5s})  ' +
              ''.join(f'{v:9.3f}' for v in R[3, 2:] / R[3, 2]))


if __name__ == '__main__':
    main(sys.argv[1])

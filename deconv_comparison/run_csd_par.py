"""run_csd.py, split across processes.

MSMT-CSD is a quadratic program per voxel and takes ~30 ms of that; at 40,000
voxels per SNR arm that is the longest single step in the package. The fits are
independent, so this runs them in a process pool. Output is byte identical to
run_csd.py -- it is the same model, the same responses and the same data, only
the loop is split.

    python3 run_csd_par.py snr30 p30 4
"""
import sys
from multiprocessing import Pool

import numpy as np
from dipy.core.gradients import gradient_table
from dipy.reconst.csdeconv import (AxSymShResponse,
                                   ConstrainedSphericalDeconvModel)
from dipy.reconst.mcsd import MultiShellDeconvModel

import binio
from run_csd import LMAX, msmt_response

_G = {}


def _init(tag, tag_resp):
    bvals = binio.load('bvals').ravel()
    bvecs = binio.load('bvecs')
    dwi = binio.load('mc_dwi_' + tag)
    S = dwi.reshape(-1, dwi.shape[-1], order='F')
    S = S / S[:, bvals == 0].mean(axis=1, keepdims=True)
    shells = np.unique(bvals)
    gtab = gradient_table(bvals * 1000.0, bvecs=bvecs)
    _G['S'] = S
    _G['model'] = MultiShellDeconvModel(gtab, msmt_response(tag_resp, shells),
                                        sh_order_max=LMAX, iso=2)


def _fit(rng):
    i0, i1 = rng
    m, S = _G['model'], _G['S']
    n = i1 - i0
    C = np.empty((n, (LMAX + 1) * (LMAX + 2) // 2))
    V = np.empty((n, 3))
    for k in range(n):
        f = m.fit(S[i0 + k])
        C[k] = f.shm_coeff
        V[k] = f.volume_fractions
    return i0, C, V


def main(tag, tag_resp='p30', nproc=4):
    nproc = int(nproc)
    bvals = binio.load('bvals').ravel()
    bvecs = binio.load('bvecs')
    dwi = binio.load('mc_dwi_' + tag)
    S = dwi.reshape(-1, dwi.shape[-1], order='F')
    S = S / S[:, bvals == 0].mean(axis=1, keepdims=True)
    N = S.shape[0]
    shells = np.unique(bvals)

    # SSST-CSD is fast enough not to be worth splitting
    b_hi = shells[-1]
    keep = (bvals == 0) | (bvals == b_hi)
    gt = gradient_table(bvals[keep] * 1000.0, bvecs=bvecs[keep])
    r_sh = binio.load(f'resp_tournier_wm_{tag_resp}')[-1].ravel()
    model = ConstrainedSphericalDeconvModel(
        gt, AxSymShResponse(1.0, r_sh), sh_order_max=LMAX)
    C = np.stack([model.fit(S[i, keep]).shm_coeff for i in range(N)])
    binio.save(f'sh_csd_{tag}', C)
    print(f'sh_csd_{tag}: {N} voxels', flush=True)
    del S, C

    edges = np.linspace(0, N, nproc * 4 + 1).astype(int)
    jobs = list(zip(edges[:-1], edges[1:]))
    Call = np.empty((N, (LMAX + 1) * (LMAX + 2) // 2))
    Vall = np.empty((N, 3))
    with Pool(nproc, initializer=_init, initargs=(tag, tag_resp)) as p:
        for i0, c, v in p.imap_unordered(_fit, jobs):
            Call[i0:i0 + len(c)] = c
            Vall[i0:i0 + len(v)] = v
    binio.save(f'sh_msmt_{tag}', Call)
    binio.save(f'vf_msmt_{tag}', Vall)
    print(f'sh_msmt_{tag}: {N} voxels', flush=True)


if __name__ == '__main__':
    main(*sys.argv[1:])

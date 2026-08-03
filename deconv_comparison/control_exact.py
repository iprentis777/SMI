"""Control: CSD and MSMT-CSD given the EXACT response instead of an estimated one.

A comparison in which MSMT-CSD comes out worst on crossings is a claim that
needs a control, because it is equally consistent with "MSMT-CSD was set up
wrong here". This runs both CSD variants a second time with the exact zonal
response of the kernel that generated the data -- an idealised delta response,
which is the best any response could possibly be -- on a subset of the
realisations. If the estimated-response numbers are a setup error, this will not
fix them; if they are the cost of estimating a response from dispersed white
matter, it will.

Writes sh_csdX_<tag> and sh_msmtX_<tag> for the first N voxels of each
condition, and a matching cond id vector, so score.py can be pointed at them.
"""
import sys
import numpy as np
from dipy.core.gradients import gradient_table
from dipy.reconst.csdeconv import (AxSymShResponse,
                                   ConstrainedSphericalDeconvModel)
from dipy.reconst.mcsd import MultiShellDeconvModel, MultiShellResponse

import binio
from kernel import response_sh

LMAX = 6
SH_CONST = 0.5 / np.sqrt(np.pi)
K_WM = [0.60, 2.0, 2.0, 0.50, 0.02]      # the kernel gen_montecarlo.m used
D_FW = 3.0
D_GM = 0.8


def main(tag, nsub=500):
    nsub = int(nsub)
    bvals = binio.load('bvals').ravel()
    bvecs = binio.load('bvecs')
    dwi = binio.load('mc_dwi_' + tag)
    S = dwi.reshape(-1, dwi.shape[-1], order='F')
    S = S / S[:, bvals == 0].mean(axis=1, keepdims=True)
    cond = binio.load('mc_cond_id_' + tag).ravel().astype(int)
    shells = np.unique(bvals)

    sub = np.concatenate([np.where(cond == c)[0][:nsub]
                          for c in np.unique(cond)])
    binio.save('mc_cond_id_' + tag + 'X', cond[sub])
    for nm in ('mc_gt_axes_', 'mc_angles_', 'mc_sh_gt6_', 'Y_smi_'):
        binio.save(nm + tag + 'X', binio.load(nm + tag))

    b_hi = shells[-1]
    keep = (bvals == 0) | (bvals == b_hi)
    gt = gradient_table(bvals[keep] * 1000.0, bvecs=bvecs[keep])
    model = ConstrainedSphericalDeconvModel(
        gt, AxSymShResponse(1.0, response_sh(b_hi, K_WM, LMAX)),
        sh_order_max=LMAX)
    C = np.stack([model.fit(S[i, keep]).shm_coeff for i in sub])
    binio.save(f'sh_csdX_{tag}X', C)
    print(f'sh_csdX_{tag}X: {len(sub)} voxels')

    R = np.zeros((len(shells), 2 + len(range(0, LMAX + 1, 2))))
    for i, b in enumerate(shells):
        R[i, 0] = np.exp(-b * D_FW) / SH_CONST
        R[i, 1] = np.exp(-b * D_GM) / SH_CONST
        R[i, 2:] = response_sh(b, K_WM, LMAX)
    resp = MultiShellResponse(R, LMAX, shells * 1000.0, S0=np.ones(len(shells)))
    mm = MultiShellDeconvModel(gradient_table(bvals * 1000.0, bvecs=bvecs),
                               resp, sh_order_max=LMAX, iso=2)
    C = np.stack([mm.fit(S[i]).shm_coeff for i in sub])
    binio.save(f'sh_msmtX_{tag}X', C)
    print(f'sh_msmtX_{tag}X: {len(sub)} voxels')


if __name__ == '__main__':
    main(sys.argv[1], *sys.argv[2:])

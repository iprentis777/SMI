"""CSD and MSMT-CSD on the compartment-split conditions, with BOTH an exact
ground-truth response and a response estimated the way the real tools estimate
it -- from the most anisotropic / most single-fibre-like voxels of a brain that
contains edema.

Variants written per tag:
  csd_b3_amp_<tag>    CSD, exact delta-fODF response from the fw=0 kernel
  csdT_b3_amp_<tag>   CSD, dwi2response-tournier response from the edema brain
  msmt_amp_<tag>      MSMT-CSD, exact per-shell responses
  msmtE_amp_<tag>     MSMT-CSD, responses estimated from the edema brain
"""
import sys
import numpy as np
from dipy.core.gradients import gradient_table
from dipy.core.sphere import Sphere
from dipy.reconst.csdeconv import (AxSymShResponse,
                                   ConstrainedSphericalDeconvModel)
from dipy.reconst.mcsd import (MultiShellDeconvModel, MultiShellResponse,
                               multi_shell_fiber_response)

import binio
from kernel import response_sh

SH_CONST = 0.5 / np.sqrt(np.pi)
LMAX = 6
D_FW = 3.0
D_GM = 0.8
TISSUE = [0.60, 2.0, 2.0, 0.50, 0.0]
BRAIN = 'edema30'


def exact_msmt_response(shell_b):
    ls = np.arange(0, LMAX + 1, 2)
    R = np.zeros((len(shell_b), 2 + len(ls)))
    for i, b in enumerate(shell_b):
        R[i, 0] = np.exp(-b * D_FW) / SH_CONST
        R[i, 1] = np.exp(-b * D_GM) / SH_CONST
        R[i, 2:] = response_sh(b, TISSUE, LMAX)
    return R


def main(tag, only=None):
    bvals = binio.load('bvals').ravel()
    bvecs = binio.load('bvecs')
    ev = binio.load('eval_dirs')
    sphere = Sphere(xyz=ev)
    dwi = binio.load('dwi_' + tag)
    S = dwi.reshape(-1, dwi.shape[-1], order='F')
    NVOX = S.shape[0]
    b0 = S[:, bvals == 0].mean(axis=1, keepdims=True)
    Sn = S / b0

    gtab_all = gradient_table(bvals * 1000.0, bvecs=bvecs)
    keep = (bvals == 0) | (bvals == 3.0)
    gt3 = gradient_table(bvals[keep] * 1000.0, bvecs=bvecs[keep])

    # ---------------------------------------------------------------- CSD
    for name, r_sh in [('csd_b3', response_sh(3.0, TISSUE, LMAX)),
                       ('csdE_b3', binio.load(f'resp_aniso_{BRAIN}')[3, 2:].ravel())]:
        if only and name not in only:
            continue
        model = ConstrainedSphericalDeconvModel(
            gt3, AxSymShResponse(1.0, np.asarray(r_sh, float)), sh_order_max=LMAX)
        amp = np.empty((ev.shape[0], NVOX))
        for i in range(NVOX):
            amp[:, i] = model.fit(Sn[i, keep]).odf(sphere)
        binio.save(f'{name}_amp_{tag}', amp)
        print(f'{name}_amp_{tag} done')

    # ----------------------------------------------------------- MSMT-CSD
    shell_b = np.unique(bvals)
    variants = {'msmt': MultiShellResponse(exact_msmt_response(shell_b), LMAX,
                                           shell_b * 1000.0,
                                           S0=np.array([1.0, 1.0, 1.0]))}
    # estimated from the most anisotropic voxels of a brain containing edema,
    # per shell and per l (no tensor step -- what MRtrix stores)
    variants['msmtE'] = MultiShellResponse(
        binio.load(f'resp_aniso_{BRAIN}'), LMAX, shell_b * 1000.0,
        S0=np.array([1.0, 1.0, 1.0]))

    for name, resp in variants.items():
        if only and name not in only:
            continue
        model = MultiShellDeconvModel(gtab_all, resp, sh_order_max=LMAX, iso=2)
        amp = np.empty((ev.shape[0], NVOX))
        vf = np.empty((NVOX, 3))
        for i in range(NVOX):
            f = model.fit(Sn[i])
            amp[:, i] = f.odf(sphere)
            vf[i] = f.volume_fractions
        binio.save(f'{name}_amp_{tag}', amp)
        binio.save(f'{name}_vf_{tag}', vf)
        print(f'{name}_amp_{tag} done')


if __name__ == '__main__':
    main(sys.argv[1], set(sys.argv[2:]) or None)

"""CSD and MSMT-CSD (dipy) on the same synthetic data SMI was fitted to.

Both are given response functions derived ANALYTICALLY from the ground-truth
kernel, so neither carries any response-estimation error. That is deliberately
generous to them: SMI has to estimate its kernel from the noisy data, they do
not. Any SMI advantage that survives this is real.

The WM response is the fw=0 tissue kernel -- what `dwi2response tournier`
would return from healthy single-fibre WM. Neither CSD nor MSMT is told that
the voxel it is being handed contains 40% free water.
"""
import sys
import numpy as np
from dipy.core.gradients import gradient_table
from dipy.core.sphere import Sphere
from dipy.reconst.csdeconv import (AxSymShResponse,
                                   ConstrainedSphericalDeconvModel)
from dipy.reconst.mcsd import MultiShellDeconvModel, MultiShellResponse

import binio
from kernel import Kell, response_sh

SH_CONST = 0.5 / np.sqrt(np.pi)          # Y_00
LMAX = 6
D_FW = 3.0
D_GM = 0.8                                # nuisance compartment; GT has no GM
TISSUE = [0.60, 2.0, 2.0, 0.50, 0.0]      # fw=0 WM response kernel


def build_msmt_response(shell_b):
    """(n_shells, 2 + LMAX//2 + 1) in multi_shell_fiber_response's layout:
    columns [csf_l0, gm_l0, wm_l0, wm_l2, ...]. Entries are SH coefficients,
    so an isotropic signal of amplitude c contributes c / Y_00."""
    ls = np.arange(0, LMAX + 1, 2)
    R = np.zeros((len(shell_b), 2 + len(ls)))
    for i, b in enumerate(shell_b):
        R[i, 0] = np.exp(-b * D_FW) / SH_CONST
        R[i, 1] = np.exp(-b * D_GM) / SH_CONST
        R[i, 2:] = response_sh(b, TISSUE, LMAX)
    return R


def main(tag):
    bvals = binio.load('bvals').ravel()            # ms/um^2
    bvecs = binio.load('bvecs')
    ev = binio.load('eval_dirs')
    sphere = Sphere(xyz=ev)

    if tag == 'clean':
        S = binio.load('S_clean_snr30')            # [NVOX x Ndwi], noise free
    else:
        dwi = binio.load('dwi_' + tag)
        S = dwi.reshape(-1, dwi.shape[-1], order='F')
    NVOX = S.shape[0]

    # normalise each voxel by its own measured b0, as every pipeline does
    b0 = S[:, bvals == 0].mean(axis=1, keepdims=True)
    Sn = S / b0

    gtab_all = gradient_table(bvals * 1000.0, bvecs=bvecs)

    # ------------------------------------------------------------ CSD
    for bshell in (3.0, 1.0):
        keep = (bvals == 0) | (bvals == bshell)
        gt = gradient_table(bvals[keep] * 1000.0, bvecs=bvecs[keep])
        resp = AxSymShResponse(1.0, response_sh(bshell, TISSUE, LMAX))
        model = ConstrainedSphericalDeconvModel(gt, resp, sh_order_max=LMAX)
        amp = np.empty((ev.shape[0], NVOX))
        for i in range(NVOX):
            amp[:, i] = model.fit(Sn[i, keep]).odf(sphere)
        name = f'csd_b{int(bshell)}_amp_{tag}'
        binio.save(name, amp)
        print(f'{name}: [{amp.shape[0]} x {amp.shape[1]}]')

    # ------------------------------------------------------- MSMT-CSD
    shell_b = np.unique(bvals)
    R = build_msmt_response(shell_b)
    resp = MultiShellResponse(R, LMAX, shell_b * 1000.0,
                              S0=np.array([1.0, 1.0, 1.0]))
    model = MultiShellDeconvModel(gtab_all, resp, sh_order_max=LMAX, iso=2)
    amp = np.empty((ev.shape[0], NVOX))
    vf = np.empty((NVOX, 3))
    for i in range(NVOX):
        f = model.fit(Sn[i])
        amp[:, i] = f.odf(sphere)
        vf[i] = f.volume_fractions
    binio.save(f'msmt_amp_{tag}', amp)
    binio.save(f'msmt_vf_{tag}', vf)
    print(f'msmt_amp_{tag}: [{amp.shape[0]} x {amp.shape[1]}]')


if __name__ == '__main__':
    main(sys.argv[1])

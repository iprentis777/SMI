"""SSST-CSD and MSMT-CSD on the Monte Carlo data, through dipy.

Both read the same `mc_dwi_<tag>` bytes that SMI was fitted to, both use the
response functions estimated by dhollander.py from the phantom, and both write
spherical harmonic coefficients in dipy's descoteaux07 basis so that peaks.py
can score all four methods with one code path.

  csd_<tag>    SSST-CSD (Tournier et al., 2007) on the b = 3 shell alone,
               with the `dwi2response tournier` response -- the pairing the
               single-shell pipeline actually uses
  msmt_<tag>   MSMT-CSD (Jeurissen et al., 2014) on all three shells, with the
               3-tissue `dwi2response dhollander` responses

Why the highest shell for SSST-CSD: single-shell CSD is normally run on the
highest available b, which is where the extra-axonal and free water signal has
decayed most and the deconvolution is best conditioned. Giving it b = 1 instead
would be a straw man.
"""
import sys
import time
import numpy as np
from dipy.core.gradients import gradient_table
from dipy.reconst.csdeconv import (AxSymShResponse,
                                   ConstrainedSphericalDeconvModel)
from dipy.reconst.mcsd import MultiShellDeconvModel, MultiShellResponse

import binio

LMAX = 6


def msmt_response(tag_resp, shells, algo='dhollander'):
    """Assemble MRtrix-style 3-tissue responses into dipy's MultiShellResponse.

    dipy expects [Nshell x (2 + Ncoef_wm)] with the two isotropic columns first
    (CSF then GM), which is the layout MultiShellResponse documents.
    """
    wm = binio.load(f'resp_{algo}_wm_{tag_resp}')
    gm = binio.load(f'resp_{algo}_gm_{tag_resp}')
    csf = binio.load(f'resp_{algo}_csf_{tag_resp}')
    R = np.zeros((len(shells), 2 + wm.shape[1]))
    R[:, 0] = csf[:, 0]
    R[:, 1] = gm[:, 0]
    R[:, 2:] = wm
    return MultiShellResponse(R, LMAX, np.asarray(shells) * 1000.0,
                              S0=np.ones(len(shells)))


def main(tag, tag_resp='p30', which=('csd', 'msmt')):
    bvals = binio.load('bvals').ravel()
    bvecs = binio.load('bvecs')
    dwi = binio.load('mc_dwi_' + tag)
    S = dwi.reshape(-1, dwi.shape[-1], order='F')
    S = S / S[:, bvals == 0].mean(axis=1, keepdims=True)
    N = S.shape[0]
    shells = np.unique(bvals)

    if 'csd' in which:
        b_hi = shells[-1]
        keep = (bvals == 0) | (bvals == b_hi)
        gt = gradient_table(bvals[keep] * 1000.0, bvecs=bvecs[keep])
        r_sh = binio.load(f'resp_tournier_wm_{tag_resp}')[-1].ravel()
        model = ConstrainedSphericalDeconvModel(
            gt, AxSymShResponse(1.0, r_sh), sh_order_max=LMAX)
        t0 = time.time()
        C = np.empty((N, (LMAX + 1) * (LMAX + 2) // 2))
        for i in range(N):
            C[i] = model.fit(S[i, keep]).shm_coeff
        binio.save(f'sh_csd_{tag}', C)
        print(f'sh_csd_{tag}: {N} voxels in {time.time()-t0:.0f} s')

    if 'msmt' in which:
        gtab = gradient_table(bvals * 1000.0, bvecs=bvecs)
        model = MultiShellDeconvModel(gtab, msmt_response(tag_resp, shells),
                                      sh_order_max=LMAX, iso=2)
        t0 = time.time()
        C = np.empty((N, (LMAX + 1) * (LMAX + 2) // 2))
        VF = np.empty((N, 3))
        for i in range(N):
            f = model.fit(S[i])
            C[i] = f.shm_coeff
            VF[i] = f.volume_fractions
        binio.save(f'sh_msmt_{tag}', C)
        binio.save(f'vf_msmt_{tag}', VF)
        print(f'sh_msmt_{tag}: {N} voxels in {time.time()-t0:.0f} s')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'p30',
         tuple(sys.argv[3:]) or ('csd', 'msmt'))

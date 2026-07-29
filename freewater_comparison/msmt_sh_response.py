"""An MRtrix-equivalent MSMT response: per shell, per l, estimated directly.

dipy's `response_from_mask_msmt` + `multi_shell_fiber_response` route the WM
response through a per-shell DIFFUSION TENSOR: three eigenvalues per shell,
from which an axially symmetric Gaussian signal is synthesised and re-projected
onto SH. That constrains the response to a Gaussian profile.

MRtrix's `dwi2response dhollander` does not do this -- it stores per-shell,
per-l response coefficients estimated by rotating each selected voxel onto a
common axis and averaging. This module does the same, so the tensor step can
be isolated as a cause rather than left as a confound.

Selection is the same tournier voxel set used for the CSD response, so the two
estimated responses come from exactly the same tissue.
"""
import numpy as np
from dipy.core.geometry import cart2sphere, vec2vec_rotmat
from dipy.core.gradients import gradient_table
from dipy.core.sphere import HemiSphere
from dipy.data import get_sphere
from dipy.direction.peaks import peaks_from_model
from dipy.reconst.csdeconv import (AxSymShResponse,
                                   ConstrainedSphericalDeconvModel)
from dipy.reconst.shm import real_sh_descoteaux_from_index
import dipy.reconst.dti as dti

import binio
from kernel import response_sh
from tournier import tournier_response

SH_CONST = 0.5 / np.sqrt(np.pi)
LMAX = 6
BRAIN = 'edema30'


def main():
    bvals = binio.load('bvals').ravel()
    bvecs = binio.load('bvecs')
    brain = binio.load(f'brain_{BRAIN}')
    shape = brain.shape[:3]
    b0 = brain[..., bvals == 0].mean(axis=-1, keepdims=True)
    data = (brain / b0).reshape(-1, len(bvals), order='F')

    # --- the same tournier voxel selection the CSD response was built from
    keep3 = (bvals == 0) | (bvals == 3.0)
    gt3 = gradient_table(bvals[keep3] * 1000.0, bvecs=bvecs[keep3])
    init = response_sh(3.0, [0.60, 2.0, 2.0, 0.50, 0.0], LMAX)
    _, idx = tournier_response(gt3, data[:, keep3], init, peak_thr=0.2)

    # --- fibre direction of each selected voxel, from the b=3 CSD peak
    r_csd = binio.load(f'resp_csd_tournier_{BRAIN}').ravel()
    model = ConstrainedSphericalDeconvModel(
        gt3, AxSymShResponse(1.0, r_csd), sh_order_max=LMAX)
    sph = HemiSphere.from_sphere(get_sphere(name='symmetric724'))
    pk = peaks_from_model(model=model, data=data[idx][:, None, None, keep3],
                          sphere=sph, relative_peak_threshold=0.2,
                          min_separation_angle=25, parallel=False)
    dirs = pk.peak_dirs.reshape(len(idx), -1, 3)[:, 0, :]

    # --- rotate-and-average, per shell
    ls = np.arange(0, LMAX + 1, 2)
    shells = np.unique(bvals)
    R = np.zeros((len(shells), 2 + len(ls)))
    for i, b in enumerate(shells):
        sel = bvals == b
        if b == 0:
            R[i, 2] = data[idx][:, sel].mean() / SH_CONST
        else:
            acc = np.zeros(len(ls))
            for j in range(len(idx)):
                Rm = vec2vec_rotmat(dirs[j], np.array([0, 0, 1.0]))
                g = (Rm @ bvecs[sel].T).T
                _, th, ph = cart2sphere(*g.T)
                B = real_sh_descoteaux_from_index(0, ls, th[:, None], ph[:, None])
                acc += np.linalg.lstsq(B, data[idx[j], sel], rcond=-1)[0]
            R[i, 2:] = acc / len(idx)

    # --- isotropic compartments from their own masks
    low = (bvals == 0) | (bvals == 1.0)
    tf = dti.TensorModel(gradient_table(bvals[low] * 1000.0,
                                        bvecs=bvecs[low])).fit(
        (brain / b0)[..., low])
    fa = tf.fa.reshape(-1, order='F')
    md = tf.md.reshape(-1, order='F')
    m_gm = (fa < 0.15) & (md < 1.5e-3)
    m_csf = (fa < 0.20) & (md > 2.0e-3)
    for i, b in enumerate(shells):
        sel = bvals == b
        R[i, 0] = data[m_csf][:, sel].mean() / SH_CONST
        R[i, 1] = data[m_gm][:, sel].mean() / SH_CONST

    binio.save('resp_msmt_sh', R)
    print('per-shell SH response (rows = b = 0/1/2/3):')
    print('           csf_l0    gm_l0    wm_l0    wm_l2    wm_l4    wm_l6')
    for i, b in enumerate(shells):
        print(f'  b={b:.0f}  ' + ''.join(f'{v:9.4f}' for v in R[i]))
    print(f'\nselected {len(idx)} voxels')
    ex = np.array([response_sh(b, [0.60, 2.0, 2.0, 0.50, 0.0], LMAX)
                   for b in shells])
    print('\nWM response / l0, estimated vs exact:')
    for i, b in enumerate(shells[1:], start=1):
        print(f'  b={b:.0f}  est {np.round(R[i,2:]/R[i,2],3)}'
              f'   exact {np.round(ex[i]/ex[i,0],3)}')


if __name__ == '__main__':
    main()

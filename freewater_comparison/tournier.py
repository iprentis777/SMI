"""The `dwi2response tournier` / dipy `recursive_response` selection, instrumented.

Algorithmically identical to dipy's `recursive_response` -- same CSD model, same
`peaks_from_model`, same second-to-first peak ratio criterion, same rotate-and-
average response accumulation -- but it returns the SELECTED VOXEL MASK as well
as the response, which dipy's version does not expose. Without the mask there is
no way to see which tissue the response was actually built from, which is the
whole question here.

dipy's default `peak_thr=0.01` (second peak below 1% of the first) selects zero
voxels on a Watson kappa=16 fODF at Lmax 6 and returns NaN; the threshold is
raised here, which is a change of operating point, not of algorithm.
"""
import numpy as np
from dipy.core.geometry import cart2sphere
from dipy.core.sphere import HemiSphere
from dipy.data import get_sphere
from dipy.direction.peaks import peaks_from_model
from dipy.reconst.csdeconv import (AxSymShResponse,
                                   ConstrainedSphericalDeconvModel)
from dipy.reconst.shm import real_sh_descoteaux_from_index
from dipy.core.geometry import vec2vec_rotmat


def tournier_response(gtab, data2d, init_response, sh_order_max=6,
                      peak_thr=0.2, n_iter=5, convergence=1e-3, n_select=300):
    """data2d is [Nvox x Ndwi]. Returns (r_sh, selected_index_array)."""
    ls = np.arange(0, sh_order_max + 1, 2)
    sphere = HemiSphere.from_sphere(get_sphere(name='symmetric724'))
    where_dwi = ~gtab.b0s_mask

    res_obj = AxSymShResponse(1.0, np.asarray(init_response, dtype=float))
    idx = np.arange(data2d.shape[0])
    resp_p = np.asarray(init_response, dtype=float)

    for _ in range(n_iter):
        model = ConstrainedSphericalDeconvModel(gtab, res_obj,
                                                sh_order_max=sh_order_max)
        pk = peaks_from_model(model=model, data=data2d[idx][:, None, None, :],
                              sphere=sphere,
                              relative_peak_threshold=peak_thr,
                              min_separation_angle=25, parallel=False)
        vals = pk.peak_values.reshape(len(idx), -1)
        dirs = pk.peak_dirs.reshape(len(idx), -1, 3)
        with np.errstate(invalid='ignore', divide='ignore'):
            ratio = np.where(vals[:, 0] > 0, vals[:, 1] / vals[:, 0], 1.0)
        # dwi2response tournier keeps a FIXED NUMBER of the most single-fibre
        # -like voxels (300 by default), not everything under a threshold.
        # Thresholding alone let 73% of the volume in and blunted the response.
        order = np.argsort(ratio)[:n_select]
        if len(order) < 10:
            break
        idx = idx[order]
        dirs = dirs[order]

        r_sh = np.zeros(len(ls))
        for j in range(len(idx)):
            R = vec2vec_rotmat(dirs[j, 0], np.array([0, 0, 1.0]))
            g = (R @ gtab.gradients.T).T[where_dwi]
            _, th, ph = cart2sphere(*g.T)
            B = real_sh_descoteaux_from_index(0, ls, th[:, None], ph[:, None])
            r_sh += np.linalg.lstsq(B, data2d[idx[j], where_dwi], rcond=-1)[0]
        r_sh /= len(idx)
        res_obj = AxSymShResponse(data2d[idx][:, gtab.b0s_mask].mean(), r_sh)
        if np.all(np.abs((resp_p - r_sh) / resp_p) < convergence):
            break
        resp_p = r_sh
    return r_sh, idx

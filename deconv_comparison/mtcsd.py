"""Multi-tissue constrained spherical deconvolution, WM + one isotropic tissue.

Needed only by the Dhollander et al. (2019) single-fibre WM selection step,
which runs `dwi2fod msmt_csd` with a WM and a CSF response and nothing else.
dipy's MultiShellDeconvModel refuses fewer than two isotropic compartments
(dipy/reconst/mcsd.py:224), so that particular fit has to be done here.

The headline MSMT-CSD results in this package do NOT come through this file --
they come from dipy's own MultiShellDeconvModel, unmodified. This exists so the
response ESTIMATION matches MRtrix, not to replace the method being compared.

Convolution convention, verified against SMI's forward model in
check_conventions.py:

    s_lm(shell) = c_lm * r_l(shell) * sqrt(4*pi/(2l+1))

with c_lm the FOD coefficients and r_l the zonal response coefficients, both in
the orthonormal real SH basis. An isotropic compartment is an FOD with only
c_00 non-zero, so its signal contribution is c_00 * r_0(shell), flat over
directions.
"""
import numpy as np

from peaks import sh_basis
from dipy.reconst.shm import sph_harm_ind_list


def design(bvals, bvecs, r_wm, r_iso, lmax, shells):
    """[Ndwi x (Ncoef + Niso)] design matrix of the multi-tissue forward model.

    r_wm   [Nshell x (lmax/2+1)] zonal WM response
    r_iso  list of [Nshell x 1] isotropic responses
    """
    _, l = sph_harm_ind_list(lmax)
    B = sh_basis(bvecs, lmax)
    A = np.zeros((len(bvals), B.shape[1] + len(r_iso)))
    for i, b in enumerate(shells):
        s = bvals == b
        rl = np.zeros(B.shape[1])
        for j, ll in enumerate(range(0, lmax + 1, 2)):
            rl[l == ll] = r_wm[i, j] * np.sqrt(4 * np.pi / (2 * ll + 1))
        A[np.ix_(s, np.arange(B.shape[1]))] = B[s] * rl
        for k, ri in enumerate(r_iso):
            A[s, B.shape[1] + k] = ri[i, 0]
    return A


class TwoTissueCSD:
    """min ||A x - S||^2 s.t. the WM FOD >= 0 on a sphere and x_iso >= 0."""

    def __init__(self, bvals, bvecs, r_wm, r_iso, lmax, reg_dirs):
        self.shells = np.unique(bvals)
        self.lmax = lmax
        self.ncoef = (lmax + 1) * (lmax + 2) // 2
        self.A = design(bvals, bvecs, r_wm, r_iso, lmax, self.shells)
        self.Breg = sh_basis(reg_dirs, lmax)
        self.niso = len(r_iso)
        self._setup()

    def _setup(self):
        import cvxpy as cp
        n = self.A.shape[1]
        self.x = cp.Variable(n)
        self.S = cp.Parameter(self.A.shape[0])
        cons = [self.Breg @ self.x[:self.ncoef] >= 0,
                self.x[self.ncoef:] >= 0]
        self.prob = cp.Problem(
            cp.Minimize(cp.sum_squares(self.A @ self.x - self.S)), cons)

    def fit(self, signal):
        self.S.value = np.asarray(signal, dtype=float)
        try:
            self.prob.solve(warm_start=True)
        except Exception:
            return np.linalg.lstsq(self.A, signal, rcond=None)[0]
        if self.x.value is None:
            return np.linalg.lstsq(self.A, signal, rcond=None)[0]
        return np.asarray(self.x.value, dtype=float)

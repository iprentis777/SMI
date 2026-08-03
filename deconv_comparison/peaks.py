"""One peak finder, applied to every method.

Every fODF in this package -- SMI constrained, SMI unconstrained, SSST-CSD,
MSMT-CSD and the band limited ground truth -- is reduced to the same
representation first (spherical harmonic coefficients in dipy's descoteaux07
basis) and then run through this file. No result can therefore come from a
different peak finder, a different sphere, or a different threshold.

Two decisions worth knowing about:

1. PEAKS ARE FOUND ON THE ANISOTROPIC PART, l >= 2. SMI's fODF is normalised to
   unit mass and carries a fixed isotropic floor of 1/(4*pi); a CSD fODF's l=0
   term is its apparent fibre density and varies. A relative peak threshold
   would therefore mean a different thing for each method. Zeroing l=0 leaves
   the positions of all local maxima unchanged and makes the threshold mean the
   same thing everywhere. Reported amplitudes are of the full fODF.

2. VERTEX QUANTISATION IS REMOVED. The coarse search runs on a 2890 vertex
   sphere (3.9 degree spacing), then each peak is refined against the continuous
   SH expansion by fitting a local quadratic in the tangent plane. Without this
   the reported angular errors would have a ~2 degree floor that has nothing to
   do with the methods.
"""
import numpy as np
from dipy.core.geometry import cart2sphere
from dipy.data import get_sphere
from dipy.direction.peaks import peak_directions
from dipy.reconst.shm import real_sh_descoteaux, sph_harm_ind_list

SPHERE = get_sphere(name='repulsion724').subdivide(n=1)
REL_THRESH = 0.30
MIN_SEP = 20.0


def sh_basis(dirs, lmax):
    """[Ndir x Ncoef] real_sh_descoteaux basis, dipy's default CSD basis."""
    _, th, ph = cart2sphere(*np.asarray(dirs, dtype=float).T)
    B, _, _ = real_sh_descoteaux(lmax, th, ph)
    return B


def l_of_coeff(lmax):
    _, l = sph_harm_ind_list(lmax)
    return l


def amplitudes(coef, B):
    """[Nvox x Ndir] amplitudes from [Nvox x Ncoef] coefficients."""
    return np.asarray(coef) @ B.T


def _tangent(d):
    """Two unit vectors spanning the tangent plane at d."""
    t = np.zeros_like(d)
    t[..., 0] = 1.0
    alt = np.abs(d[..., 0]) > 0.9
    t[alt] = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(d, t)
    e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)
    e2 = np.cross(d, e1)
    return e1, e2


def refine(coef, dirs, lmax, steps=(4.0, 1.0, 0.25)):
    """Refine peak directions against the continuous SH expansion.

    coef  [N x Ncoef], dirs [N x 3] -- one starting direction per row. Returns
    (refined dirs, amplitude at them). A six point stencil in the tangent plane
    determines a full 2-D quadratic exactly; its stationary point is the next
    iterate, rejected if it leaves the stencil (which happens on saddle-like
    neighbourhoods, where the coarse vertex is already the best estimate).
    """
    coef = np.atleast_2d(np.asarray(coef, dtype=float))
    d = np.asarray(dirs, dtype=float).copy()
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    for h_deg in steps:
        h = np.radians(h_deg)
        e1, e2 = _tangent(d)
        s = 1.0 / np.sqrt(2.0)
        offs = [(0, 0), (h, 0), (-h, 0), (0, h), (0, -h), (h * s, h * s)]
        P = np.stack([d + a * e1 + b * e2 for a, b in offs], axis=1)   # [N,6,3]
        P /= np.linalg.norm(P, axis=2, keepdims=True)
        B = sh_basis(P.reshape(-1, 3), lmax).reshape(len(d), 6, -1)
        f = np.einsum('nkc,nc->nk', B, coef)
        # design of the quadratic a + b1 u + b2 v + c11 u^2/2 + c22 v^2/2 + c12 uv
        M = np.array([[1, u, v, 0.5 * u * u, 0.5 * v * v, u * v]
                      for u, v in offs], dtype=float)
        q = np.linalg.solve(M, f.T).T                                  # [N,6]
        b = q[:, 1:3]
        Hm = np.stack([np.stack([q[:, 3], q[:, 5]], -1),
                       np.stack([q[:, 5], q[:, 4]], -1)], -2)          # [N,2,2]
        det = Hm[:, 0, 0] * Hm[:, 1, 1] - Hm[:, 0, 1] * Hm[:, 1, 0]
        ok = np.abs(det) > 1e-12
        step = np.zeros_like(b)
        if ok.any():
            step[ok] = -np.linalg.solve(Hm[ok], b[ok][..., None])[..., 0]
        ok &= np.all(np.abs(step) <= h, axis=1)
        nd = d + step[:, :1] * e1 + step[:, 1:2] * e2
        nd /= np.linalg.norm(nd, axis=1, keepdims=True)
        d = np.where(ok[:, None], nd, d)
    B = sh_basis(d, lmax)
    amp = np.einsum('nc,nc->n', B, coef)
    return d, amp


def find_peaks(coef, lmax, B_sphere=None, rel_thresh=REL_THRESH,
               min_sep=MIN_SEP, max_peaks=5, chunk=2000):
    """Peaks of every fODF in `coef`.

    Returns a list of (dirs [k x 3], amps [k]) with k <= max_peaks, sorted by
    amplitude. Search is on the anisotropic part; the returned amplitude is of
    the full fODF, so the isotropic term is included.
    """
    coef = np.atleast_2d(np.asarray(coef, dtype=float))
    if B_sphere is None:
        B_sphere = sh_basis(SPHERE.vertices, lmax)
    aniso = coef.copy()
    aniso[:, 0] = 0.0

    out = []
    for i0 in range(0, len(coef), chunk):
        blk = aniso[i0:i0 + chunk]
        A = blk @ B_sphere.T
        for j in range(len(blk)):
            a = A[j]
            if not np.isfinite(a).all() or a.max() <= 0:
                out.append((np.zeros((0, 3)), np.zeros(0)))
                continue
            d, v, _ = peak_directions(a, SPHERE,
                                      relative_peak_threshold=rel_thresh,
                                      min_separation_angle=min_sep)
            if len(d) == 0:
                out.append((np.zeros((0, 3)), np.zeros(0)))
                continue
            out.append((d[:max_peaks], v[:max_peaks]))

    # refine every peak of every voxel in one batch
    counts = np.array([len(d) for d, _ in out])
    if counts.sum() == 0:
        return out
    rows = np.repeat(np.arange(len(out)), counts)
    starts = np.concatenate([d for d, _ in out if len(d)], axis=0)
    rd, ramp = refine(coef[rows], starts, lmax)
    res, k = [], 0
    for n in counts:
        order = np.argsort(ramp[k:k + n])[::-1]
        res.append((rd[k:k + n][order], ramp[k:k + n][order]))
        k += n
    return res


def angle(u, v):
    """Angle in degrees between undirected axes."""
    return np.degrees(np.arccos(np.clip(np.abs(np.dot(u, v)), 0, 1)))


def match_error(pk_dirs, true_axes):
    """Greedy match of found peaks to true axes.

    Returns (n_peaks, mean angular error over matched axes, matched_all).
    Peaks are visited in amplitude order, so the largest peak is matched first.
    """
    if len(pk_dirs) == 0:
        return 0, np.nan, False
    used, errs = set(), []
    for a in true_axes:
        best, bi = np.inf, None
        for i, p in enumerate(pk_dirs):
            if i in used:
                continue
            e = angle(p, a)
            if e < best:
                best, bi = e, i
        if bi is None:
            return len(pk_dirs), np.nan, False
        used.add(bi)
        errs.append(best)
    return len(pk_dirs), float(np.mean(errs)), True


def acc(u, v):
    """Angular correlation coefficient (Anderson 2005) over l >= 2.

    Scale free and basis independent, so it compares fODFs whose amplitudes
    are not on the same scale -- which SMI's and CSD's are not.
    """
    u = np.asarray(u)[..., 1:]
    v = np.asarray(v)[..., 1:]
    num = (u * v).sum(-1)
    den = np.sqrt((u * u).sum(-1) * (v * v).sum(-1))
    return np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)

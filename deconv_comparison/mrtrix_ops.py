"""The MRtrix3 primitives that `dwi2response dhollander` is built out of.

Each function here reproduces one MRtrix command closely enough that the
selection it drives is the same selection. Where the reproduction is not exact
the docstring says so, because the whole point of implementing the algorithm
rather than approximating it is that the response functions handed to CSD and
MSMT-CSD are the ones those tools would actually be given.

  optimal_threshold  `mrthreshold` with no options (Ridgway et al., 2009)
  threshold_top      `mrthreshold -top N` / `-bottom N`
  amp2response       `amp2response`, including the non-negativity and
                     monotonicity constraints
"""
import numpy as np


def optimal_threshold(x):
    """`mrthreshold` default: the threshold maximising image/mask correlation.

    MRtrix minimises

        cost(t) = -cov(x, x > t) / (std(x) * std(x > t))

    over t (core/filter/optimal_threshold.h). MRtrix searches with a golden
    section search between min and max, which finds a local optimum; this
    evaluates the cost at every candidate threshold instead, which is exact and
    cheap at these array sizes. Returns t such that `x > t` is the selection.
    """
    x = np.asarray(x, dtype=float).ravel()
    x = x[np.isfinite(x)]
    n = x.size
    if n == 0:
        return np.nan
    xs = np.sort(x)
    total = xs.sum()
    # candidate thresholds sit between consecutive distinct values; selecting
    # `x > t` with t = xs[i] keeps the n-i-1 largest
    csum = np.cumsum(xs[::-1])[::-1]          # csum[i] = sum of xs[i:]
    k = n - np.arange(n)                      # count of xs[i:]
    # for threshold just below xs[i], the mask is xs[i:]
    mean_x = total / n
    std_x = np.sqrt(max((xs @ xs) / n - mean_x**2, 0.0))
    if std_x <= 0:
        return float(xs[0])
    mean_xy = csum / n
    frac = k / n
    cov = mean_xy - frac * mean_x
    var_m = frac - frac**2                    # variance of the 0/1 mask
    with np.errstate(invalid='ignore', divide='ignore'):
        corr = np.where(var_m > 0, cov / (std_x * np.sqrt(var_m)), -np.inf)
    corr[0] = -np.inf                         # the all-ones mask has zero var
    i = int(np.argmax(corr))
    # a threshold strictly below xs[i] and at or above xs[i-1]
    return float(np.nextafter(xs[i], -np.inf))


def threshold_top(x, n, bottom=False, ignorezero=True):
    """`mrthreshold -top N` (or `-bottom N`): boolean mask of the N extremes.

    `-ignorezero` excludes exact zeros from the candidate set, which is how
    dhollander.py uses it -- the images it thresholds are zero outside the mask
    of interest.
    """
    x = np.asarray(x, dtype=float).ravel()
    ok = np.isfinite(x)
    if ignorezero:
        ok &= x != 0
    idx = np.where(ok)[0]
    if idx.size == 0:
        return np.zeros(x.shape, bool)
    n = int(min(n, idx.size))
    order = np.argsort(x[idx], kind='stable')
    pick = idx[order[:n]] if bottom else idx[order[::-1][:n]]
    m = np.zeros(x.shape, bool)
    m[pick] = True
    return m


def _zonal_basis(theta, lmax):
    """[len(theta) x (lmax/2+1)] matrix of Y_l0(theta), even l only."""
    x = np.cos(np.asarray(theta, dtype=float))
    P = {0: np.ones_like(x),
         2: 1.5 * x**2 - 0.5,
         4: (35 * x**4 - 30 * x**2 + 3) / 8,
         6: (231 * x**6 - 315 * x**4 + 105 * x**2 - 5) / 16,
         8: (6435 * x**8 - 12012 * x**6 + 6930 * x**4 - 1260 * x**2 + 35) / 128}
    ls = np.arange(0, lmax + 1, 2)
    return np.stack([np.sqrt((2 * l + 1) / (4 * np.pi)) * P[l] for l in ls], axis=1)


def amp2response(amps, gdirs, fibre_dirs, lmax, isotropic=False,
                 constrained=True):
    """`amp2response`: zonal harmonic response from a set of selected voxels.

    amps        [Nvox x Ndir] signal amplitudes of one shell
    gdirs       [Ndir x 3] gradient directions of that shell
    fibre_dirs  [Nvox x 3] the fibre direction of each voxel (MRtrix uses the
                DTI principal eigenvector here, not an FOD peak)
    lmax        maximum even order
    isotropic   fit l = 0 only, as `amp2response -isotropic` does

    Every voxel is rotated onto a common axis and all of them enter ONE least
    squares problem, which is what amp2response does -- not a per-voxel fit
    followed by an average.

    `constrained` adds amp2response's two constraints, on the same grid it uses
    (1 degree steps from 0 to 90, cmd/amp2response.cpp:343-354): the response
    must be non-negative, and its amplitude must INCREASE from the fibre
    direction out to the orthogonal plane. The direction of that second
    constraint is the whole point -- the diffusion signal is smallest along the
    fibre -- and getting its sign wrong collapses the fit to an isotropic
    response, which is how it was caught. Without cvxpy the fit falls back to
    plain least squares.
    """
    amps = np.atleast_2d(np.asarray(amps, dtype=float))
    if isotropic:
        return np.array([amps.mean() * np.sqrt(4 * np.pi)])

    fibre_dirs = np.asarray(fibre_dirs, dtype=float)
    gdirs = np.asarray(gdirs, dtype=float)
    # undirected: the response is antipodally symmetric, so |cos| is the angle
    cos_t = np.clip(np.abs(fibre_dirs @ gdirs.T), 0.0, 1.0)
    theta = np.arccos(cos_t)                       # [Nvox x Ndir]
    B = _zonal_basis(theta.ravel(), lmax)          # [Nvox*Ndir x Ncoef]
    y = amps.ravel()

    if not constrained:
        return np.linalg.lstsq(B, y, rcond=None)[0]

    try:
        import cvxpy as cp
    except ImportError:
        return np.linalg.lstsq(B, y, rcond=None)[0]

    th_g = np.radians(np.arange(91.0))
    G = _zonal_basis(th_g, lmax)
    c = cp.Variable(B.shape[1])
    cons = [G @ c >= 0, cp.diff(G @ c) >= 0]       # non-negative, non-decreasing
    prob = cp.Problem(cp.Minimize(cp.sum_squares(B @ c - y)), cons)
    try:
        prob.solve()
    except Exception:
        return np.linalg.lstsq(B, y, rcond=None)[0]
    if c.value is None:
        return np.linalg.lstsq(B, y, rcond=None)[0]
    return np.asarray(c.value, dtype=float)

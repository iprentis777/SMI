"""Shared loading and peak extraction.

Peak extraction is deliberately ONE implementation applied to all methods, on
one shared sphere, so no result can come from a different peak finder.
"""
import numpy as np
from dipy.data import get_sphere
from dipy.direction.peaks import peak_directions

import binio

SPHERE = get_sphere(name='repulsion724').subdivide(n=2)
REL_THRESH = 0.35
MIN_SEP = 25.0


def cond_names(suf=''):
    with open(binio.DATA + f'/cond_names{suf}.txt') as f:
        return [l.strip() for l in f if l.strip()]


def gt_axes(suf=''):
    """[3 x 3 x NCOND] true fibre axes, NaN-padded."""
    return binio.load('gt_axes' + suf)


def cond_id(suf=''):
    return binio.load('cond_id' + suf).ravel().astype(int)


def peaks_of(amp_col):
    """Peaks of one voxel's ODF.

    The relative-amplitude threshold is applied to the ANISOTROPIC part
    (ODF minus its spherical mean). SMI's fODF carries a fixed isotropic floor
    of 1/(4*pi) that CSD's does not, and without removing it a relative
    threshold would mean something different for each method. Local maxima
    locations are unchanged by the subtraction, so principal-peak directions
    are unaffected either way.
    """
    a = amp_col - amp_col.mean()
    a = np.clip(a, 0, None)
    if a.max() <= 0:
        return np.zeros((0, 3)), np.zeros(0)
    d, v, _ = peak_directions(a, SPHERE,
                              relative_peak_threshold=REL_THRESH,
                              min_separation_angle=MIN_SEP)
    return d, v


def ang(u, v):
    """Angle in degrees between undirected axes."""
    c = np.clip(abs(float(np.dot(u, v))), 0, 1)
    return np.degrees(np.arccos(c))


def principal_error(amp_col, axes):
    """Angle between the largest peak and the nearest true fibre axis."""
    d, _ = peaks_of(amp_col)
    if len(d) == 0:
        return np.nan
    return min(ang(d[0], a) for a in axes)


def matched_error(amp_col, axes):
    """Greedy match of found peaks to true axes.

    Returns (n_peaks, mean angular error over matched axes, all_matched).
    """
    d, _ = peaks_of(amp_col)
    if len(d) == 0:
        return 0, np.nan, False
    errs, used = [], set()
    for a in axes:
        best, bi = np.inf, None
        for i, p in enumerate(d):
            if i in used:
                continue
            e = ang(p, a)
            if e < best:
                best, bi = e, i
        if bi is None:
            return len(d), np.nan, False
        used.add(bi)
        errs.append(best)
    return len(d), float(np.mean(errs)), True


def peak_amp(amp_col):
    """Raw ODF amplitude at its maximum (floor included)."""
    return float(amp_col.max())

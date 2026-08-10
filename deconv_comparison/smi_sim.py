"""The SMI simulation's forward model, in Python.

This is the Python mirror of the ground truth, kernel and forward convolution
that `notebooks/smi_manuscript_60deg.m` builds in Octave. It exists so the CSD
arms can be driven from a Python notebook while provably simulating the *same*
experiment as the SMI arm, rather than a similar one.

WHAT IS MIRRORED, AND FROM WHERE

    uniform_sphere_dirs   helpers/fODF_modulation_helpers.m  uniform_sphere_dirs
    get_even_SH           SMI.m:3194                         SMI.get_even_SH
    watson_amp            helpers/fODF_modulation_helpers.m  watson_amp
    fodf_to_plm           helpers/fODF_modulation_helpers.m  fodf_to_plm
    Kell                  SMI.m:2415   RotInv_Kell_wFW_b_beta_TE_numerical
    forward_signal        helpers/fODF_modulation_helpers.m  forward_signal
    response_zh           helpers/SMI_response_helpers.m     kernel_zh
    load_protocol         mc_config.m                        load_protocol
    pick_grid             mc_config.m                        pick_grid
    rotate_about          mc_config.m                        rotate_about

Every one of those is checked against the Octave original by `check_python_vs_octave.py`,
which runs both and reports max|err|. Nothing here is trusted because it looks
right; the notebook re-runs that comparison as a CHECK on every execution.

WHAT IS DELIBERATELY *NOT* MIRRORED

The noise realisations. The .m file draws from Octave's legacy generator
(`randn('seed',...)`) and this module draws from numpy's PCG64. Those streams
cannot be made bit-identical without reimplementing one inside the other, and
pretending otherwise would be the kind of assumption this package exists to
avoid. What is guaranteed instead:

  * the NOISE-FREE signal is identical to machine precision, which is the claim
    that both sides simulate the same experiment;
  * the noise *process* is the same -- Rician, sigma = 1/SNR against S0 = 1 --
    and is verified by recovering sigma from the data on both sides;
  * the noisy signal this module generates is EXPORTED, so the SMI arm can be
    run on the very same realisations if common random numbers are wanted.

CONVENTIONS. b is in ms/um^2 throughout, diffusivities in um^2/ms. fODF
coefficients are in SMI's normalised convention p_00 = 1, in which the fODF
integrates to 1 in every voxel. `plm` covers l = 2..Lmax and excludes the l = 0
term, matching `out.plm`; `sh` vectors include l = 0 and are what MRtrix stores.
"""
import os
import numpy as np
from scipy.special import lpmv

ROOT = os.path.dirname(os.path.abspath(__file__))

# The 200-point Gauss-Legendre rule SMI uses, mapped to [0,1]. SMI.m carries the
# nodes and weights as a hardcoded table; numpy generates the same rule, and
# check_python_vs_octave.py measures the difference rather than assuming it away.
_GL_N = 200


def _gauss_legendre_01(n=_GL_N):
    """Nodes and weights of the n-point Gauss-Legendre rule on [0,1].

    Weights sum to 1, matching the table at SMI.m:2470.
    """
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


# Legendre polynomials in x, even orders only, matching SMI.m:2543-2553 and
# helpers/SMI_response_helpers.m legendre_P. Written out rather than recursed so
# they are visibly the same expressions as the Octave side.
_P_EVEN = {
    0: lambda x: np.ones_like(x),
    2: lambda x: 1.5 * x**2 - 0.5,
    4: lambda x: (35 * x**4 - 30 * x**2 + 3) / 8,
    6: lambda x: (231 * x**6 - 315 * x**4 + 105 * x**2 - 5) / 16,
    8: lambda x: (6435 * x**8 - 12012 * x**6 + 6930 * x**4 - 1260 * x**2 + 35) / 128,
}

LMAX_KERNEL_CEILING = 8      # SMI's K_l are undefined above l = 8


# ---------------------------------------------------------------- geometry
def uniform_sphere_dirs(N):
    """N nearly uniform directions on the sphere, as an [N x 3] array.

    Golden angle spiral, usable as an equal-area quadrature grid. Mirrors
    helpers/fODF_modulation_helpers.m uniform_sphere_dirs.
    """
    ga = np.pi * (3 - np.sqrt(5))
    k = np.arange(N, dtype=float)
    z = 1 - 2 * (k + 0.5) / N
    r = np.sqrt(np.maximum(0, 1 - z**2))
    th = ga * k
    return np.column_stack([r * np.cos(th), r * np.sin(th), z])


def rotate_about(n, deg):
    """A unit vector `deg` degrees from n, in an arbitrary but FIXED plane.

    Mirrors mc_config.m rotate_about, so every arm that asks for the same angle
    gets the same pair of axes.
    """
    n = np.asarray(n, dtype=float).ravel()
    n = n / np.linalg.norm(n)
    t = np.array([0.0, 0.0, 1.0])
    if abs(n @ t) > 0.9:
        t = np.array([1.0, 0.0, 0.0])
    e = t - (t @ n) * n
    e = e / np.linalg.norm(e)
    m = np.cos(np.deg2rad(deg)) * n + np.sin(np.deg2rad(deg)) * e
    return m / np.linalg.norm(m)


def pick_grid(N):
    """A 3D grid holding exactly N voxels with every dimension > 1.

    Mirrors mc_config.m pick_grid. SMI.vectorize takes a different branch if any
    spatial dimension is a singleton, so a simulation laid out as [N 1 1] is
    silently a different code path. The most cube-like factorisation wins.
    """
    divisors = [d for d in range(1, N + 1) if N % d == 0]
    best = None
    for a in divisors:
        if not (1 < a < N):
            continue
        m = N // a
        for b in [d for d in range(1, m + 1) if m % d == 0]:
            if not (1 < b < m):
                continue
            c = m // b
            if c > 1:
                cand = tuple(sorted([a, b, c]))
                if best is None or (max(cand) - min(cand)) < (max(best) - min(best)):
                    best = cand
    if best is None:
        raise ValueError(
            'pick_grid: %d does not factor into three integers > 1. '
            'Choose NREP so that NREP * NCOND does.' % N)
    return list(best)


# ------------------------------------------------------------ SH machinery
def sh_degrees(Lmax):
    """The l of every even-order SH coefficient, in SMI's ordering."""
    return np.concatenate([np.full(2 * l + 1, l) for l in range(0, Lmax + 1, 2)])


def sh_orders(Lmax):
    """The m of every even-order SH coefficient, in SMI's ordering."""
    return np.concatenate([np.arange(-l, l + 1) for l in range(0, Lmax + 1, 2)])


def n_coef(Lmax):
    return (Lmax // 2 + 1) * (Lmax + 1)


def get_even_SH(dirs, Lmax, CS_phase=0):
    """[N x ncoef] real even-order SH basis, mirroring SMI.get_even_SH (SMI.m:3194).

    CS_phase selects SMI's two conventions. At CS_phase = 0 this basis IS
    MRtrix's, which is what lets an SMI fODF and an MRtrix FOD be scored by the
    same lines of code; at 1 the two differ by (-1)^m, a 180 degree rotation
    about z of every fODF.

    Note the phase bookkeeping, which is easy to get backwards: MATLAB's and
    Octave's `legendre` already include the Condon-Shortley (-1)^m, and so does
    scipy's `lpmv`. SMI's CS_phase = 1 branch multiplies by a further (-1)^m,
    which CANCELS it. So the flag does not mean what its name suggests, and the
    only safe thing to do is reproduce the arithmetic exactly, which is what is
    below.
    """
    dirs = np.asarray(dirs, dtype=float)
    if dirs.shape[1] != 3:
        dirs = dirs.T
    x, y, z = dirs[:, 0], dirs[:, 1], dirs[:, 2]
    PHI = np.arctan2(y, x)                      # cart2sph azimuth
    THETA = np.arccos(np.clip(z, -1.0, 1.0))    # pi/2 - elevation, for unit dirs

    l_all = sh_degrees(Lmax)
    m_all = sh_orders(Lmax)
    am = np.abs(m_all)

    from scipy.special import gammaln
    # sqrt((2l+1)/(4pi) * (l-|m|)!/(l+|m|)!), via log-gamma so l = 8, m = 8 does
    # not overflow the factorial the way a literal transcription would.
    K_lm = np.sqrt((2 * l_all + 1) / (4 * np.pi)
                   * np.exp(gammaln(l_all - am + 1) - gammaln(l_all + am + 1)))

    extra = np.ones_like(K_lm)
    extra[m_all != 0] = np.sqrt(2)
    if CS_phase:
        extra = extra * (-1.0) ** m_all

    ct = np.cos(THETA)
    Y = np.empty((dirs.shape[0], len(l_all)))
    for i, (l, m) in enumerate(zip(l_all, m_all)):
        P = lpmv(abs(m), l, ct)                 # includes (-1)^|m|, as MATLAB does
        if m > 0:
            phi = np.cos(m * PHI)
        elif m == 0:
            phi = np.ones_like(PHI)
        else:
            phi = np.sin(-m * PHI)
        Y[:, i] = extra[i] * K_lm[i] * phi * P
    return Y


# ------------------------------------------------------------- the fODF
def watson_amp(dirs_q, n, kappa):
    """Unnormalised Watson amplitudes about axis n. kappa <= 0 gives isotropic.

    Mirrors helpers/fODF_modulation_helpers.m watson_amp.
    """
    dirs_q = np.asarray(dirs_q, dtype=float)
    if kappa <= 0:
        return np.ones(dirs_q.shape[0])
    n = np.asarray(n, dtype=float).ravel()
    return np.exp(kappa * (dirs_q @ (n / np.linalg.norm(n))) ** 2)


def fodf_to_plm(fodf, dirs_q, Lmax, CS_phase=0):
    """Project a sampled fODF onto normalised plm, in SMI's p_00 = 1 convention.

    Returns l = 2..Lmax only, matching `out.plm`. Mirrors
    helpers/fODF_modulation_helpers.m fodf_to_plm.
    """
    fodf = np.asarray(fodf, dtype=float).ravel()
    Y = get_even_SH(dirs_q, Lmax, CS_phase)
    wq = 4 * np.pi / dirs_q.shape[0]            # equal-area quadrature weight
    f = fodf / (fodf.sum() * wq)                # so that \int fODF dOmega = 1
    flm = (Y.T @ f) * wq
    L = sh_degrees(Lmax)
    plm_all = flm / np.sqrt((2 * L + 1) / (4 * np.pi))
    return plm_all[1:] / plm_all[0]


def plm_to_sh(plm, Lmax):
    """SH coefficients (l = 0 included) from a normalised plm vector.

        f_lm = plm * sqrt((2l+1)/(4pi)),   p_00 = 1

    This is the vector MRtrix stores, and the one every glyph is drawn from.
    """
    L = sh_degrees(Lmax)
    full = np.concatenate([[1.0], np.asarray(plm, dtype=float).ravel()])
    return full * np.sqrt((2 * L + 1) / (4 * np.pi))


def watson_zonal_pl(kappa, Lmax, ndir=None, dirs_q=None, CS_phase=0):
    """Normalised zonal coefficients p_l (p_0 = 1) of a Watson about +z.

    These are what turns a delta response into a DISPERSION-MATCHED one: a
    single fibre population in this simulation is not a delta, it is a Watson at
    this kappa, so the response of one fibre is the kernel convolved with it.

    Computed by the same sampled-grid projection the ground truth uses, so the
    response and the truth cannot disagree about what a Watson is.
    `watson_zonal_pl_exact` computes the same numbers by 1-D quadrature with no
    spherical harmonics anywhere, and the two are compared in the notebook.
    """
    if dirs_q is None:
        dirs_q = uniform_sphere_dirs(ndir if ndir else 3000)
    f = watson_amp(dirs_q, [0, 0, 1], kappa)
    plm = fodf_to_plm(f, dirs_q, Lmax, CS_phase)
    sh = np.concatenate([[1.0], plm])
    L = sh_degrees(Lmax)
    M = sh_orders(Lmax)
    return sh[M == 0]                            # the m = 0 entries, l = 0,2,..


def watson_zonal_pl_exact(kappa, Lmax, n=2000):
    """The same p_l as watson_zonal_pl, by 1-D Gauss-Legendre and no harmonics.

    For an axially symmetric fODF f(x), x = cos(theta), the normalised zonal
    coefficients in the p_00 = 1 convention are

        p_l = int_0^1 f(x) P_l(x) dx / int_0^1 f(x) dx

    (the [-1,0] half adds nothing new: f is even in x and so is every even P_l).
    An independent route to the same numbers, so a mistake in the SH projection
    cannot hide inside the response as well as inside the truth.
    """
    x, w = _gauss_legendre_01(n)
    f = np.exp(kappa * x**2) if kappa > 0 else np.ones_like(x)
    denom = w @ f
    return np.array([(w @ (f * _P_EVEN[l](x))) / denom
                     for l in range(0, Lmax + 1, 2)])


# ------------------------------------------------------------- the kernel
def Kell(ell, b, kernel, D_FW=3.0, n=_GL_N):
    """K_l(b) for kernel = [f, Da, Depar, Deperp, fw], b in ms/um^2.

    Mirrors SMI.RotInv_Kell_wFW_b_beta_TE_numerical (SMI.m:2415) for LTE
    (beta = 1) and a single TE, which is what this whole package simulates:

      K_l(b) = int_0^1 [ f e^{-b Da x^2}
                       + (1-f-fw) e^{-b(Deperp + (Depar-Deperp) x^2)}
                       + fw e^{-b D_FW} ] P_l(x) dx

    There is no (2l+1) prefactor. The free water term is carried at every l and
    removes itself for even l >= 2, because int_0^1 P_l dx = 0 there.
    """
    f, Da, Depar, Deperp, fw = np.asarray(kernel, dtype=float).ravel()[:5]
    b = np.atleast_1d(np.asarray(b, dtype=float))
    x, w = _gauss_legendre_01(n)
    P = _P_EVEN[ell](x)
    dpar = Depar - Deperp
    integ = (f * np.exp(-np.outer(b, Da * x**2))
             + (1 - f - fw) * np.exp(-np.outer(b, Deperp + dpar * x**2))
             + fw * np.exp(-np.outer(b, D_FW * np.ones_like(x))))
    out = (integ * (w * P)).sum(axis=1)
    if ell > 0:
        out = np.where(b < 1e-6, 0.0, out)      # SMI.m, the b ~ 0 special case
    return out


def Kell_matrix(kernel, b, Lmax, D_FW=3.0):
    """[Nb x (Lmax/2+1)] rotational invariants, one column per even l.

    Mirrors helpers/SMI_response_helpers.m kernel_Kell.
    """
    if Lmax > LMAX_KERNEL_CEILING:
        raise ValueError('SMI kernel invariants are undefined above l = %d'
                         % LMAX_KERNEL_CEILING)
    b = np.atleast_1d(np.asarray(b, dtype=float))
    return np.column_stack([Kell(l, b, kernel, D_FW)
                            for l in range(0, Lmax + 1, 2)])


def response_zh(kernel, b, Lmax, D_FW=3.0, pl=None, S0=1.0):
    """Zonal harmonic response, [Nb x (Lmax/2+1)], in MRtrix's response format.

    With `pl` omitted this is the DELTA response, the one SMI's own convention
    implies for a single fibre:

        r_l(b) = K_l(b) sqrt((2l+1) 4 pi)

    `pl` supplies the normalised zonal coefficients of the fibre population the
    response is meant to describe, giving a DISPERSION-MATCHED response:

        r_l(b) = K_l(b) p_l sqrt((2l+1) 4 pi)

    Pass `watson_zonal_pl(kappa, Lmax)` to get the response of one Watson
    population at the same kappa the ground truth uses. That is the honest
    single-fibre response for this simulation: a "fibre" here is a dispersed
    population, not a delta, and a response estimated from real white matter
    would have absorbed that dispersion too.
    """
    K = Kell_matrix(kernel, b, Lmax, D_FW)
    l = np.arange(0, Lmax + 1, 2)
    R = S0 * K * np.sqrt((2 * l + 1) * 4 * np.pi)[None, :]
    if pl is not None:
        R = R * np.asarray(pl, dtype=float).ravel()[None, :]
    return R


def zh_profile(r, theta):
    """Amplitude of an axially symmetric function with zonal coefficients r."""
    r = np.asarray(r, dtype=float).ravel()
    x = np.cos(theta)
    A = np.zeros_like(np.asarray(theta, dtype=float))
    for i, l in enumerate(range(0, 2 * len(r), 2)):
        A = A + r[i] * np.sqrt((2 * l + 1) / (4 * np.pi)) * _P_EVEN[l](x)
    return A


# ------------------------------------------------------- forward convolution
def forward_signal(plm, kernel, bvals, bvecs, Lmax, CS_phase=0, D_FW=3.0):
    """Noise-free S(u)/S0 for one voxel, as a [Ndwi] array.

        S(u)/S0 = sum_lm K_l(b) p_lm Y_lm(u) sqrt((2l+1) 4 pi)

    Exactly the construction SMI.get_plm_from_S_and_kernel inverts. Mirrors
    helpers/fODF_modulation_helpers.m forward_signal.
    """
    bvals = np.asarray(bvals, dtype=float).ravel()
    L = sh_degrees(Lmax)
    N_l = np.sqrt((2 * L + 1) * 4 * np.pi)
    Y = get_even_SH(bvecs, Lmax, CS_phase)
    K = Kell_matrix(kernel, bvals, Lmax, D_FW)      # [Ndwi x (Lmax/2+1)]
    Kmat = K[:, L // 2]                              # expand to one column per lm
    coef = np.concatenate([[1.0], np.asarray(plm, dtype=float).ravel()])
    return (Kmat * (Y * N_l)) @ coef


def forward_signal_direct(fodf_q, dirs_q, kernel, bvals, bvecs, D_FW=3.0):
    """The same signal with no spherical harmonics anywhere.

        S(u) = sum_q w_q fODF(n_q) K(u . n_q),   sum_q w_q = 1

    Direct numerical convolution over the quadrature grid. This is the check
    that catches a basis or normalisation error, which a self-consistent SH
    round trip would hide. It is NOT band limited, so it does not agree with
    forward_signal to machine precision -- the residual IS the ground truth's
    band-limiting error.
    """
    f, Da, Depar, Deperp, fw = np.asarray(kernel, dtype=float).ravel()[:5]
    bvals = np.asarray(bvals, dtype=float).ravel()
    cosang = np.asarray(bvecs, dtype=float) @ np.asarray(dirs_q, dtype=float).T
    Kmat = (f * np.exp(-bvals[:, None] * (Da * cosang**2))
            + (1 - f - fw) * np.exp(-bvals[:, None] * (Deperp + (Depar - Deperp) * cosang**2))
            + fw * np.exp(-bvals[:, None] * D_FW))
    w = np.asarray(fodf_q, dtype=float).ravel()
    w = w / w.sum()
    return Kmat @ w


# ------------------------------------------------------------------ noise
def rician(S, sigma, rng):
    """Rician magnitude noise: |S + sigma n1 + i sigma n2|.

    At sigma = 0 this returns sqrt(S^2) = S, i.e. the noise-free signal exactly.
    The draws still happen so the stream has the same shape at every SNR, which
    matches what the .m file does.
    """
    S = np.asarray(S, dtype=float)
    n1 = rng.standard_normal(S.shape)
    n2 = rng.standard_normal(S.shape)
    return np.sqrt((S + sigma * n1) ** 2 + (sigma * n2) ** 2)


# --------------------------------------------------------------- the protocol
def load_protocol(fname='hcp_real_3shell.txt', b0_snap=0.0, verbose=True):
    """Read protocol/<fname>. Returns (bvals in ms/um^2, bvecs [N x 3] unit).

    Mirrors mc_config.m load_protocol, including its warning, because the
    warning is the point: a real .bvec is written at seven significant figures
    and is unit only to about 1e-6. Anything that reads g(3) as cos(theta)
    without normalising inherits that, and at Lmax 8 it degrades the
    zonal-response identity from ~1e-15 to ~5e-7.

    `b0_snap` treats any shell below that b as exactly 0, which restores the
    exact identity S(0)/S0 = 1 -- a real .bval records EFFECTIVE b, so its
    "b = 0" volumes are b = 5 s/mm^2.
    """
    path = fname if os.path.isabs(fname) else os.path.join(ROOT, 'protocol', fname)
    rows = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith('%') or ln.startswith('#'):
                continue
            rows.append([float(v) for v in ln.split()])
    arr = np.asarray(rows, dtype=float)
    bvals, bvecs = arr[:, 0].copy(), arr[:, 1:4].copy()

    nrm = np.sqrt((bvecs**2).sum(axis=1))
    e_raw = np.abs(nrm - 1).max()
    if e_raw > 1e-12 and verbose:
        bar = '!' * 72
        print('\n%s' % bar)
        print('WARNING  %s: gradient directions are not unit vectors.' % os.path.basename(path))
        print('         max | |g| - 1 | = %.3e over %d volumes.' % (e_raw, len(bvals)))
        print('         They are being normalised. Left uncorrected this breaks any')
        print('         calculation that reads g(3) as cos(theta) -- measured at Lmax 8,')
        print('         the zonal response vs forward model identity degrades from')
        print('         ~1e-15 to ~5e-7. Check the source of the .bvec if this is large.')
        print('%s\n' % bar)
    bvecs = bvecs / nrm[:, None]

    if b0_snap > 0:
        n_snap = int(((bvals > 0) & (bvals < b0_snap)).sum())
        if n_snap and verbose:
            print('   note: %d volumes with 0 < b < %g snapped to exactly b = 0'
                  % (n_snap, b0_snap))
            print('         (acquired as b = %g ms/um^2 = %.0f s/mm^2, effective b'
                  % (bvals[bvals < b0_snap].max(), 1000 * bvals[bvals < b0_snap].max()))
            print('          from the imaging gradients)')
        bvals[bvals < b0_snap] = 0.0
    return bvals, bvecs


def group_shells(bvals, tol=0.1):
    """Bin b values into shells the way a human would, ascending.

    Returns (shell_b, shell_id) with shell_id 0-based into shell_b. The real
    protocol's b values jitter -- 18 distinct values across 4 shells -- so this
    has to cluster rather than match exact values. The notebook checks the
    result against SMI's own binning and against MRtrix's.
    """
    bvals = np.asarray(bvals, dtype=float).ravel()
    order = np.argsort(bvals)
    shell_id = np.full(bvals.shape, -1, dtype=int)
    centres = []
    for i in order:
        placed = False
        for k, c in enumerate(centres):
            if abs(bvals[i] - c) <= tol * max(c, 1.0):
                shell_id[i] = k
                placed = True
                break
        if not placed:
            centres.append(bvals[i])
            shell_id[i] = len(centres) - 1
    centres = np.array([bvals[shell_id == k].mean() for k in range(len(centres))])
    remap = np.argsort(centres)
    inv = np.empty_like(remap)
    inv[remap] = np.arange(len(remap))
    return centres[remap], inv[shell_id]


# ---------------------------------------------------------- MRtrix image IO
def write_mrtrix(name, A, grad=None, datatype='Float32LE', vox=None):
    """Write an MRtrix .mih header plus .dat raw data.

    The Python mirror of mrtrix_io.m write_image, and the only MRtrix behaviour
    this package implements locally. Data is written with strides 1,2,3,4, i.e.
    plain column-major with axis 0 fastest, and `mrinfo` is used in the notebook
    to check that MRtrix reads back what we think we wrote.
    """
    A = np.asarray(A)
    sz = list(A.shape) + [1] * max(0, 3 - A.ndim)
    if vox is None:
        vox = [1] * len(sz)
    base, _ = os.path.splitext(name)
    dtypes = {'Float32LE': '<f4', 'Float64LE': '<f8', 'UInt8': 'u1'}
    if datatype not in dtypes:
        raise ValueError('write_mrtrix: unsupported datatype %s' % datatype)

    with open(base + '.dat', 'wb') as fh:
        fh.write(np.asfortranarray(A.reshape(sz, order='F')).astype(dtypes[datatype]).tobytes(order='F'))

    with open(base + '.mih', 'w') as fh:
        fh.write('mrtrix image\n')
        fh.write('dim: %s\n' % ','.join(str(int(s)) for s in sz))
        fh.write('vox: %s\n' % ','.join('%g' % v for v in vox))
        fh.write('layout: %s\n' % ','.join('+%d' % k for k in range(len(sz))))
        fh.write('datatype: %s\n' % datatype)
        fh.write('transform: 1,0,0,0\n')
        fh.write('transform: 0,1,0,0\n')
        fh.write('transform: 0,0,1,0\n')
        if grad is not None:
            for row in np.asarray(grad, dtype=float):
                fh.write('dw_scheme: %.10g,%.10g,%.10g,%.10g\n' % tuple(row))
        fh.write('file: %s 0\n' % os.path.basename(base + '.dat'))
        fh.write('END\n')
    return base + '.mih'


def read_mrtrix(name):
    """Read a .mif or .mih image as an ndarray indexed [x, y, z, ...].

    Honours whatever `layout` the file declares, so an MRtrix output written
    with any strides comes back in the same index order the Octave side uses.
    A thin re-export of the reader in mrtrix_io.py, kept here so a notebook
    imports one module rather than two.
    """
    import mrtrix_io
    return mrtrix_io.read(name)


def write_response(path, R):
    """Write [Nshell x Ncoef] zonal coefficients as an MRtrix response .txt."""
    R = np.atleast_2d(np.asarray(R, dtype=float))
    with open(path, 'w') as fh:
        for row in R:
            fh.write(' '.join('%.8g' % v for v in row) + '\n')
    return path


def read_response(path):
    """Read an MRtrix response .txt, skipping '#' comment lines."""
    rows = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            rows.append([float(v) for v in ln.split()])
    return np.asarray(rows, dtype=float)


# ------------------------------------------------- Octave/Python binary exchange
def _binio_dir(data_dir=None):
    d = data_dir or os.path.join(ROOT, 'data')
    os.makedirs(d, exist_ok=True)
    return d


def bin_save(name, arr, data_dir=None):
    """Write a float64 array where binio.m can read it, with a shape sidecar."""
    d = _binio_dir(data_dir)
    arr = np.asarray(arr, dtype=np.float64)
    with open(os.path.join(d, name + '.bin'), 'wb') as fh:
        fh.write(np.asfortranarray(arr).tobytes(order='F'))
    with open(os.path.join(d, name + '.shape'), 'w') as fh:
        fh.write(' '.join(str(s) for s in arr.shape))


def bin_load(name, data_dir=None):
    """Read an array written by bin_save or by binio.m."""
    d = _binio_dir(data_dir)
    with open(os.path.join(d, name + '.shape')) as fh:
        shape = tuple(int(s) for s in fh.read().split())
    a = np.fromfile(os.path.join(d, name + '.bin'), dtype=np.float64)
    return a.reshape(shape, order='F')

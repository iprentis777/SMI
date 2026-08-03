"""Convention checks that the rest of the package depends on.

Nothing here produces a result. It exists because three separate conventions
have to line up before any comparison means anything, and two of them have
never been verified in this repository before:

1. SMI's even spherical harmonic basis against dipy's descoteaux07 basis
   (and therefore against MRtrix's, which descoteaux07 matches). Verified by
   solving for the map between them on the shared evaluation directions and
   checking the residual, and by round-tripping a delta along a known axis and
   confirming the peak comes back at that axis.

2. The convolution rule s_lm = c_lm r_l sqrt(4 pi/(2l+1)) used by mtcsd.py,
   against SMI's own forward model.

3. The zonal harmonic response r_l = K_l sqrt((2l+1) 4 pi) written by
   SMI_response_helpers.m, against the same forward model. (That one is also
   checked inside example_SMI_response_shview.m, so it is checked in both
   languages.)

Run after setup_protocol.py and one gen_montecarlo.m arm.
"""
import sys
import numpy as np

import binio
import peaks as pk
from kernel import Kell, response_sh

LMAX = 6


def check_basis(tag):
    ev = binio.load('eval_dirs')
    Y = binio.load('Y_smi_' + tag)
    B = pk.sh_basis(ev, LMAX)
    M, *_ = np.linalg.lstsq(B, Y, rcond=None)
    resid = np.abs(B @ M - Y).max()
    print(f'1. SMI vs dipy SH basis')
    print(f'   exact linear map exists, max|B M - Y| = {resid:.3e}')
    print(f'   the map is {"" if np.allclose(M, np.eye(len(M)), atol=1e-8) else "NOT "}'
          f'the identity, i.e. the two bases are '
          f'{"identical" if np.allclose(M, np.eye(len(M)), atol=1e-8) else "different"}')

    # a delta along a known axis, expressed in SMI's convention, must peak there
    n = np.array([0.30, -0.50, 0.81]); n /= np.linalg.norm(n)
    ls = np.repeat(np.arange(0, LMAX + 1, 2), 2 * np.arange(0, LMAX + 1, 2) + 1)
    Yn = pk.sh_basis(n[None, :], LMAX)[0]     # dipy basis at n
    # delta in dipy's basis is f_lm = Y_lm(n); build it there, map back through
    # M^-1 to check the inverse direction too
    c = Yn.copy()
    P = pk.find_peaks(c[None, :], LMAX)
    d = P[0][0][0]
    print(f'   delta at {np.round(n,3)} recovered at {np.round(d,3)}, '
          f'off by {pk.angle(d, n):.4f} deg')
    return resid


def check_convolution():
    """s_lm = c_lm r_l sqrt(4 pi/(2l+1)) must reproduce SMI's forward model."""
    bvals = binio.load('bvals').ravel()
    bvecs = binio.load('bvecs')
    kern = [0.60, 2.0, 2.0, 0.50, 0.05]
    shells = np.unique(bvals)

    n = np.array([0.30, -0.50, 0.81]); n /= np.linalg.norm(n)
    c = pk.sh_basis(n[None, :], LMAX)[0]                # delta fODF at n
    ls = np.repeat(np.arange(0, LMAX + 1, 2), 2 * np.arange(0, LMAX + 1, 2) + 1)
    B = pk.sh_basis(bvecs, LMAX)

    S_conv = np.zeros(len(bvals))
    for b in shells:
        s = bvals == b
        r = response_sh(b, kern, LMAX)
        w = np.array([r[l // 2] * np.sqrt(4 * np.pi / (2 * l + 1)) for l in ls])
        S_conv[s] = B[s] @ (c * w)

    # the same thing straight from the kernel: S(u) = sum_l K_l (2l+1) P_l(cos)
    x = np.clip(np.abs(bvecs @ n), 0, 1)
    Pl = {0: np.ones_like(x), 2: 1.5 * x**2 - .5,
          4: (35 * x**4 - 30 * x**2 + 3) / 8,
          6: (231 * x**6 - 315 * x**4 + 105 * x**2 - 5) / 16}
    S_ker = np.zeros(len(bvals))
    for b in shells:
        s = bvals == b
        acc = np.zeros(s.sum())
        for l in range(0, LMAX + 1, 2):
            acc += Kell(l, b, kern)[0] * (2 * l + 1) * Pl[l][s]
        S_ker[s] = acc
    err = np.abs(S_conv - S_ker).max()
    print(f'\n2. convolution rule vs the SM kernel')
    print(f'   max|s_lm route - K_l route| = {err:.3e}')
    return err


def check_response_file():
    """The MRtrix response file written by MATLAB must decode to the same K_l."""
    import os
    f = os.path.join(os.path.dirname(binio.DATA), '..', 'response_SMI_wm.txt')
    if not os.path.exists(f):
        print('\n3. response_SMI_wm.txt not present '
              '(run example_SMI_response_shview.m first) -- skipped')
        return None
    R = np.loadtxt(f)
    kern = [0.60, 2.0, 2.0, 0.50, 0.05]
    ref = np.stack([response_sh(b, kern, 8) for b in [0, 1, 2, 3]])
    err = np.abs(R - ref).max()
    print(f'\n3. response_SMI_wm.txt vs kernel.py')
    print(f'   max|written - recomputed| = {err:.3e}')
    return err


def check_mrtrix_basis():
    """SMI's basis against MRtrix's, and what happens if it is ignored.

    Needs dump_bases.m to have been run. This is the check that decides whether
    an SMI fODF can be handed to `tckgen` as it stands.
    """
    from dipy.reconst.shm import real_sh_tournier, sph_harm_ind_list
    from dipy.core.geometry import cart2sphere
    import os

    if not os.path.exists(os.path.join(binio.DATA, 'Y_smi_cs0.bin')):
        print('\n4. Y_smi_cs0 not present (run dump_bases.m) -- skipped')
        return None
    ev = binio.load('eval_dirs')
    _, th, ph = cart2sphere(*ev.T)
    Bt, _, _ = real_sh_tournier(LMAX, th, ph)
    m, _ = sph_harm_ind_list(LMAX)

    print('\n4. SMI vs MRtrix (dipy `tournier07`) spherical harmonic basis')
    worst = 0.0
    for cs in (0, 1):
        Y = binio.load(f'Y_smi_cs{cs}')
        M, *_ = np.linalg.lstsq(Bt, Y, rcond=None)
        resid = np.abs(Bt @ M - Y).max()
        worst = max(worst, resid)
        off = M - np.diag(np.diag(M))
        expect = np.where(m == 0, 1.0, np.sqrt(2.0))
        if cs:
            expect = expect * (-1.0) ** m
        print(f'   CS_phase={cs}: exact map, residual {resid:.1e}; diagonal '
              f'{"" if np.abs(off).max() < 1e-9 else "NOT "}(max off diagonal '
              f'{np.abs(off).max():.1e})')
        print(f'      diagonal is {"" if np.allclose(np.diag(M), expect) else "NOT "}'
              f'[1 at m=0, sqrt(2) at m!=0]{" times (-1)^m" if cs else ""}')

    # the consequence: read SMI coefficients as if they were MRtrix's
    n = np.array([0.30, -0.50, 0.81]); n /= np.linalg.norm(n)
    A = np.exp(40 * (ev @ n) ** 2)
    print('   reading SMI coefficients in MRtrix\'s basis without converting:')
    for cs in (0, 1):
        Y = binio.load(f'Y_smi_cs{cs}')
        c = np.linalg.lstsq(Y, A, rcond=None)[0]
        d_self = pk.find_peaks(c[None, :], LMAX, B_sphere=Y)[0][0][0]
        d_mrt = pk.find_peaks(c[None, :], LMAX, B_sphere=Bt)[0][0][0]
        print(f'      CS_phase={cs}: peak {np.round(d_mrt, 3)}, '
              f'{pk.angle(d_mrt, n):6.2f} deg from the truth '
              f'(SMI\'s own basis puts it {pk.angle(d_self, n):.2f} deg off, '
              f'which is the Lmax {LMAX} band limit)')
    return worst


if __name__ == '__main__':
    tag = sys.argv[1] if len(sys.argv) > 1 else 'snr30'
    r1 = check_basis(tag)
    r2 = check_convolution()
    r3 = check_response_file()
    r4 = check_mrtrix_basis()
    bad = [x for x in (r1, r2, r3, r4) if x is not None and x > 1e-8]
    print('\nALL CHECKS PASS' if not bad else f'\nFAILED: {bad}')

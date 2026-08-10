"""Measure the Python forward model against the Octave one, array by array.

This is the file that turns "the Python notebook simulates the same experiment
as smi_manuscript_60deg.m" from a claim into a number. Run dump_reference.m
first; it writes the Octave side into data/.

    octave-cli --no-gui -q dump_reference.m
    python3 check_python_vs_octave.py

Every row of the table below recomputes one array in Python and reports
max|err| against Octave's. The tolerances are not uniform, and the differences
are the interesting part:

  * The SH basis, the Watson amplitudes, the fODF projection, the kernel
    invariants and the noise-free signal must agree to ROUNDING. They are the
    same closed-form expressions evaluated in two languages, so anything above
    ~1e-12 is a real disagreement, not accumulation.

  * K_l is the one place the two sides use genuinely different arithmetic:
    SMI.m carries a hardcoded 200-point Gauss-Legendre table and numpy
    generates the rule. They should still agree to ~1e-14; if this row is the
    only one that moves, the table is the reason.

  * The harmonic-vs-direct residual is NOT expected to be zero on either side.
    It is the ground truth's band-limiting error, and what is compared is that
    the two sides get the SAME band-limiting error.

Exit status is 0 only if every row passes.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smi_sim as S                                          # noqa: E402

# The manuscript configuration, and it must match dump_reference.m exactly.
K_WM = [0.60, 2.0, 2.0, 0.50, 0.02]
D_FW, KAPPA, ANGLE = 3.0, 16.0, 60.0
AXIS1 = np.array([0.30, -0.50, 0.81])
AXIS1 = AXIS1 / np.linalg.norm(AXIS1)
LMAX_GT, CS_PHASE = 8, 0
PROTOCOL, B0_SNAP = 'hcp_real_3shell.txt', 0.05
NDIR_Q, NDIR_E = 3000, 1500

RESULTS = []


def check(name, got, want, tol, note='', rel=False):
    """Compare two arrays and record the result.

    `rel` divides the error by the scale of the reference. Use it wherever the
    quantity is not O(1): the Watson amplitudes are exp(kappa) ~ 9e6 at
    kappa = 16, and an absolute tolerance there is really a demand for 22
    significant digits.
    """
    got = np.asarray(got, dtype=float)
    want = np.asarray(want, dtype=float)
    if got.shape != want.shape:
        want = want.reshape(got.shape)
    err = float(np.max(np.abs(got - want))) if got.size else 0.0
    if rel:
        scale = float(np.max(np.abs(want))) or 1.0
        err = err / scale
    ok = err <= tol
    RESULTS.append((name, err, tol, ok, note))
    print('  %-34s max|%serr| = %9.2e   tol %7.1e   %s%s'
          % (name, 'rel ' if rel else '', err, tol,
             'ok' if ok else '** FAILED **', ('   ' + note) if note else ''))
    return ok


def main():
    data = os.path.join(S.ROOT, 'data')
    if not os.path.exists(os.path.join(data, 'ref_S_clean.shape')):
        print('data/ has no Octave reference. Run first:\n'
              '    octave-cli --no-gui -q dump_reference.m')
        return 2

    L = S.bin_load
    print('\n=== Python forward model vs the Octave original ===')
    print('config: %d deg, kappa %g, kernel %s, Lmax %d, CS_phase %d\n'
          % (ANGLE, KAPPA, K_WM, LMAX_GT, CS_PHASE))

    # ---- the protocol. Read through the same file by the same rules.
    bvals, bvecs = S.load_protocol(PROTOCOL, b0_snap=B0_SNAP, verbose=False)
    # Not exact, and the reason is worth recording: Octave's textscan and
    # Python's float() disagree by 1 ulp on decimal strings like "2.99" --
    # measured, 80 of 288 values, 1.5e-16 relative. Both are reading the same
    # text; their decimal-to-binary conversions round differently in the last
    # bit. It propagates into exp(-b D) as ~1e-16 and nothing downstream can
    # see it, but a tolerance of 0 here would fail forever for no reason.
    check('protocol b values', bvals, L('ref_bvals').ravel(), 1e-15,
          'textscan vs float(): 1 ulp', rel=True)
    check('protocol directions (unit)', bvecs, L('ref_bvecs'), 1e-15)

    # ---- the direction sets. A mismatch here would silently move every
    # quadrature sum below, so it is checked before anything uses them.
    dq = S.uniform_sphere_dirs(NDIR_Q)
    de = S.uniform_sphere_dirs(NDIR_E)
    check('quadrature directions', dq, L('ref_dq'), 1e-15)
    check('evaluation directions', de, L('ref_de'), 1e-15)

    # ---- the SH basis, which is the single most load-bearing convention here.
    check('SH basis Lmax 8 (quadrature)', S.get_even_SH(dq, LMAX_GT, CS_PHASE),
          L('ref_Y_dq_L8'), 1e-13)
    check('SH basis Lmax 6 (evaluation)', S.get_even_SH(de, 6, CS_PHASE),
          L('ref_Y_de_L6'), 1e-13)
    check('SH basis Lmax 8 (gradients)', S.get_even_SH(bvecs, LMAX_GT, CS_PHASE),
          L('ref_Y_g_L8'), 1e-13)
    check('SH basis at CS_phase = 1', S.get_even_SH(de, 6, 1),
          L('ref_Y_de_L6_cs1'), 1e-13, 'the other convention')

    # ---- the ground truth geometry and fODF.
    axis2 = S.rotate_about(AXIS1, ANGLE)
    check('second fibre axis', axis2, L('ref_axis2').ravel(), 1e-15)

    # Relative, because a Watson at kappa = 16 peaks at exp(16) ~ 8.9e6 and an
    # absolute tolerance would be asking for far more digits than a double has.
    fodf = (S.watson_amp(dq, AXIS1, KAPPA) + S.watson_amp(dq, axis2, KAPPA))
    check('Watson mixture, sampled', fodf, L('ref_fodf_q').ravel(), 1e-14,
          'peaks at exp(kappa) ~ 9e6', rel=True)

    plm_gt = S.fodf_to_plm(fodf, dq, LMAX_GT, CS_PHASE)
    check('ground truth plm', plm_gt, L('ref_plm_gt').ravel(), 1e-12)

    # ---- the Watson zonal coefficients, which make the response dispersion
    # matched. Checked against Octave AND against an independent 1-D quadrature.
    pl = S.watson_zonal_pl(KAPPA, LMAX_GT, dirs_q=dq, CS_phase=CS_PHASE)
    check('Watson zonal p_l', pl, L('ref_watson_pl').ravel(), 1e-12)
    check('Watson p_l, SH vs 1-D quadrature', pl,
          S.watson_zonal_pl_exact(KAPPA, LMAX_GT), 2e-4,
          'independent route, band limited')

    # ---- the kernel. The one row where the two sides genuinely differ in
    # method: a hardcoded quadrature table against a generated one.
    check('kernel invariants K_l(b)', S.Kell_matrix(K_WM, bvals, LMAX_GT, D_FW),
          L('ref_Kell'), 1e-13, 'table vs generated Gauss-Legendre')
    check('delta response r_l(b)', S.response_zh(K_WM, bvals, LMAX_GT, D_FW),
          L('ref_resp_delta'), 1e-12)

    # ---- the noise-free signal. This is THE claim: same experiment.
    S_clean = S.forward_signal(plm_gt, K_WM, bvals, bvecs, LMAX_GT, CS_PHASE, D_FW)
    check('noise-free signal (harmonic)', S_clean, L('ref_S_clean').ravel(), 1e-12,
          '<-- the identical-experiment claim')

    S_direct = S.forward_signal_direct(fodf, dq, K_WM, bvals, bvecs, D_FW)
    check('noise-free signal (direct)', S_direct, L('ref_S_direct').ravel(), 1e-12)

    # The band-limiting residual is a property of the truth, not an error. What
    # matters is that both sides measure the SAME one.
    e_py = float(np.max(np.abs(S_clean - S_direct)))
    e_oc = float(np.max(np.abs(L('ref_S_clean').ravel() - L('ref_S_direct').ravel())))
    check('band-limiting residual agrees', [e_py], [e_oc], 1e-12,
          'py %.2e vs oct %.2e, not zero on either side' % (e_py, e_oc))

    # ---- pick_grid, whose factorisation decides the SMI.vectorize code path.
    for N in (81, 2187, 3000, 27):
        check('pick_grid(%d)' % N, S.pick_grid(N),
              L('ref_grid_%d' % N).ravel(), 0.0)

    nfail = sum(1 for r in RESULTS if not r[3])
    print('\n%d of %d comparisons pass.' % (len(RESULTS) - nfail, len(RESULTS)))
    if nfail:
        print('FAILED rows:')
        for name, err, tol, ok, _ in RESULTS:
            if not ok:
                print('   %-34s %9.2e > %7.1e' % (name, err, tol))
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())

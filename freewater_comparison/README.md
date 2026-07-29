# Free water comparison: SMI vs CSD vs MSMT-CSD

Simulation behind `../REPORT_fODF_freewater.md`. Measures how the regularized
SMI deconvolution (no modulation) recovers a WM fODF in 40% free water,
against single-shell CSD and MSMT-CSD as implemented in dipy.

The DWI is synthesised once by SMI's own forward model and written to disk;
Octave and Python read the same bytes, and every method's fODF is evaluated on
the same sphere with the same peak finder. See section 8 of the report for the
run order and section 2 for why the design is fair.

Requires Octave (`statistics`, `image`) and Python
(`numpy scipy dipy cvxpy matplotlib`). Intermediate arrays land in `data/`,
which is not tracked.

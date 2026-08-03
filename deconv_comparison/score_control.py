"""Score the exact-response control written by control_exact.py.

Same peak finder, same ground truth, same tables as score.py -- only the method
list differs, so the numbers sit on the same scale as the main comparison.
"""
import sys

import score

if __name__ == '__main__':
    tag = sys.argv[1]
    score.METHODS = [('SSST-CSD exact resp', 'sh_csdX_%s'),
                     ('MSMT-CSD exact resp', 'sh_msmtX_%s')]
    score.summarise(tag, sys.argv[2] if len(sys.argv) > 2 else tag)

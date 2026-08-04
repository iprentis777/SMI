"""The acquisition protocol.

dipy is used here and nowhere else in the package: only to generate
electrostatically repulsed gradient directions. Everything downstream --
Octave/SMI and MRtrix -- reads the table this writes, so every method sees
identical gradient directions.

The protocol is the one Jeurissen et al. (2014) used to introduce MSMT-CSD:
three shells at b = 1, 2, 3 ms/um^2 with 90 directions each plus b = 0, i.e. an
HCP-like acquisition. b is carried in ms/um^2 throughout (b = 1 is
1000 s/mm^2); dipy is handed b*1000.
"""
import numpy as np
from dipy.core.sphere import HemiSphere, disperse_charges
from dipy.data import get_sphere

import binio

# ---------------------------------------------------------------- protocol
SHELLS = [(0.0, 18), (1.0, 90), (2.0, 90), (3.0, 90)]


def uniform_dirs(n, seed):
    """n electrostatically-repulsed directions on the hemisphere."""
    r = np.random.default_rng(seed)
    v = r.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    hs, _ = disperse_charges(HemiSphere(xyz=v), 400)
    return hs.vertices


bvals, bvecs = [], []
for i, (b, n) in enumerate(SHELLS):
    d = np.zeros((n, 3)) if b == 0 else uniform_dirs(n, 1000 + i)
    bvals.append(np.full(n, b))
    bvecs.append(d)
bvals = np.concatenate(bvals)
bvecs = np.concatenate(bvecs, axis=0)

# ------------------------------------------------------- evaluation sphere
# A dense direction set, used only by dump_bases.m and sweep_nonneg.m to
# evaluate spherical harmonics on a common grid. Peaks are NOT found on it:
# every peak in this package comes from MRtrix's sh2peaks, which does its own
# Newton search on the continuous SH expansion.
sphere = get_sphere(name='repulsion724').subdivide(n=1)
verts = sphere.vertices

binio.save('bvals', bvals)
binio.save('bvecs', bvecs)
binio.save('eval_dirs', verts)

nn = np.degrees(np.arccos(np.clip(
    np.sort(verts @ verts.T, axis=1)[:, -2], -1, 1))).mean()
print(f'protocol    : {len(bvals)} volumes, shells {SHELLS}')
print(f'eval sphere : {verts.shape[0]} vertices, mean nn spacing {nn:.2f} deg')

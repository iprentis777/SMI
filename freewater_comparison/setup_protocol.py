"""Build the acquisition protocol and the common evaluation sphere.

Everything downstream (Octave/SMI and Python/dipy) reads these files, so all
three methods see identical gradient directions and are evaluated on identical
directions. That removes protocol and peak-finder differences as confounds.
"""
import numpy as np
from dipy.core.sphere import HemiSphere, disperse_charges
from dipy.data import get_sphere

import binio

rng = np.random.default_rng(20240729)

# ---------------------------------------------------------------- protocol
# HCP-like multi-shell, b in ms/um^2 (b=1 is 1000 s/mm^2).
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
    if b == 0:
        d = np.zeros((n, 3))
    else:
        d = uniform_dirs(n, 1000 + i)
    bvals.append(np.full(n, b))
    bvecs.append(d)
bvals = np.concatenate(bvals)
bvecs = np.concatenate(bvecs, axis=0)

# ------------------------------------------------------- evaluation sphere
# repulsion724 subdivided twice: dense enough that vertex quantisation
# (~1.9 deg spacing) is negligible against the angular errors we expect.
sphere = get_sphere(name='repulsion724').subdivide(n=2)
verts = sphere.vertices

binio.save('bvals', bvals)
binio.save('bvecs', bvecs)
binio.save('eval_dirs', verts)

print(f'protocol : {len(bvals)} volumes, shells '
      f'{[(b, n) for b, n in SHELLS]}')
print(f'eval sphere : {verts.shape[0]} vertices, '
      f'mean nn spacing {np.degrees(np.arccos(np.sort(verts @ verts.T, axis=1)[:, -2])).mean():.2f} deg')

"""Flat float64 binary exchange between Octave and Python.

Octave writes column-major, so every array is read back with order='F'.
Shapes live in a sidecar '<name>.shape' text file so neither side has to
guess.
"""
import os
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ROOT, 'data')
os.makedirs(DATA, exist_ok=True)


def save(name, arr):
    arr = np.asarray(arr, dtype=np.float64)
    with open(os.path.join(DATA, name + '.bin'), 'wb') as f:
        f.write(np.asfortranarray(arr).tobytes(order='F'))
    with open(os.path.join(DATA, name + '.shape'), 'w') as f:
        f.write(' '.join(str(s) for s in arr.shape))


def load(name):
    with open(os.path.join(DATA, name + '.shape')) as f:
        shape = tuple(int(s) for s in f.read().split())
    a = np.fromfile(os.path.join(DATA, name + '.bin'), dtype=np.float64)
    return a.reshape(shape, order='F')


def exists(name):
    return os.path.exists(os.path.join(DATA, name + '.bin'))

"""Read MRtrix images (.mif / .mih) in Python.

The mirror of the reader in mrtrix_io.m, and the only way the scoring sees
anything MRtrix produced. `layout` is honoured, so an image written with any
strides comes back in the same index order as the Octave side uses.
"""
import os
import numpy as np

_DT = {'Float32LE': '<f4', 'Float32BE': '>f4', 'Float32': '<f4',
       'Float64LE': '<f8', 'Float64BE': '>f8', 'Float64': '<f8',
       'UInt8': 'u1', 'Bit': 'u1',
       'Int16LE': '<i2', 'UInt16LE': '<u2'}


def read(name):
    """Return the image as an ndarray indexed [x, y, z, ...]."""
    if not os.path.splitext(name)[1]:
        name += '.mih'
    d = os.path.dirname(name)
    dim = layout = None
    datatype = ''
    datafile = ''
    offset = 0
    with open(name, 'rb') as f:
        while True:
            raw = f.readline()
            if not raw:                       # a .mih has no END line
                break
            ln = raw.decode('utf-8', 'replace').strip()
            if ln == 'END':
                offset_after_header = f.tell()
                break
            if ':' not in ln:
                continue
            key, val = ln.split(':', 1)
            key, val = key.strip(), val.strip()
            if key == 'dim':
                dim = [int(x) for x in val.split(',')]
            elif key == 'layout':
                layout = [x.strip() for x in val.split(',')]
            elif key == 'datatype':
                datatype = val
            elif key == 'file':
                parts = val.split()
                datafile = parts[0]
                if len(parts) > 1:
                    offset = int(parts[1])
        else:
            offset_after_header = f.tell()

    if datafile == '.':
        datafile = name
        if offset == 0:
            offset = offset_after_header
    else:
        datafile = os.path.join(d, datafile)

    if datatype not in _DT:
        raise ValueError(f'mrtrix_io: unsupported datatype {datatype}')
    a = np.fromfile(datafile, dtype=np.dtype(_DT[datatype]),
                    count=int(np.prod(dim)), offset=offset).astype(np.float64)

    # `layout` gives, per axis, its position in the file ordering; the axis
    # whose entry is 0 varies fastest.
    sign = np.array([1 - 2 * (s[0] == '-') for s in layout])
    pos = np.array([int(s[1:]) for s in layout])
    file_order = np.argsort(pos)                 # axes fastest to slowest
    a = a.reshape([dim[k] for k in file_order], order='F')
    a = np.transpose(a, np.argsort(file_order))
    for k in np.where(sign < 0)[0]:
        a = np.flip(a, axis=int(k))
    return np.ascontiguousarray(a)


def read_peaks(name, npeaks=None):
    """sh2peaks output -> (directions [Nvox x P x 3], amplitudes [Nvox x P]).

    sh2peaks stores each peak as its direction scaled by its amplitude, and
    fills unfound peaks with NaN.
    """
    a = read(name)
    v = a.reshape(-1, a.shape[-1], order='F')
    P = v.shape[1] // 3
    if npeaks:
        P = min(P, npeaks)
    v = v[:, :3 * P].reshape(-1, P, 3)
    amp = np.linalg.norm(v, axis=2)
    with np.errstate(invalid='ignore', divide='ignore'):
        d = v / amp[..., None]
    return d, amp

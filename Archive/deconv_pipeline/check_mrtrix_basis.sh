#!/bin/bash
# Is SMI's spherical harmonic basis MRtrix's?
#
# Answered with MRtrix, not with a reimplementation of it. Two tests:
#
#   1. `sh2amp` is handed an SH image whose 28 voxels each hold one unit basis
#      vector, so its output IS MRtrix's basis matrix on the given directions.
#      That is compared against SMI.get_even_SH with CS_phase 0 and 1.
#   2. `sh2peaks` is handed a sharp fODF along a known axis, written in SMI's
#      convention at each CS_phase, and asked where the peak is.
#
# The second test is the one that matters for tractography: it is exactly what
# happens when out.plm is written to a .mif and passed to tckgen.
set -e
cd "$(dirname "$0")"
W=mrtrix/basis_check
mkdir -p $W

octave-cli --no-gui -q --eval "
run('$(pwd)/oct_path.m');
M = mrtrix_io(); IO = binio(); HF = fODF_modulation_helpers();
ev = IO.load('eval_dirs'); ev = ev(1:200,:);
IO.save('basischeck_dirs', ev);
E = zeros(28,1,1,28); for k=1:28, E(k,1,1,k) = 1; end
M.write('$W/idsh', E);
fid = fopen('$W/dirs.txt','w'); fprintf(fid,'%.12g %.12g %.12g\n', ev'); fclose(fid);
n = [0.30 -0.50 0.81]; n = n/norm(n);
dq = HF.dirs(3000);
for cs = [0 1]
  plm = HF.mixture_plm(HF.watson(dq,n,40), dq, 6, cs);
  L = repelem(0:2:6, 2*(0:2:6)+1)';
  c = [1; plm(:)].*sqrt((2*L+1)/(4*pi));
  V = zeros(2,2,2,28); for k=1:28, V(:,:,:,k) = c(k); end
  M.write(sprintf('$W/fod_cs%d',cs), V);
end
" 2>/dev/null | grep -v '^$' || true

sh2amp -quiet $W/idsh.mih $W/dirs.txt $W/idamp.mih -force

echo "=== 1. SMI.get_even_SH vs MRtrix's basis (via sh2amp)"
octave-cli --no-gui -q --eval "
run('$(pwd)/oct_path.m');
M = mrtrix_io(); IO = binio();
ev = IO.load('basischeck_dirs');
A  = M.read('$W/idamp.mih');
Bm = squeeze(A(:,1,1,:))';                       % [Ndir x 28], MRtrix's basis
for cs = [0 1]
  Ys = SMI.get_even_SH(ev, 6, cs);
  T  = Bm\Ys; T(abs(T) < 1e-6) = 0;
  offdiag = max(max(abs(T - diag(diag(T)))));
  m = repelem(-(0:2:6), 1); mm = [];
  for l = 0:2:6, mm = [mm -l:l]; end
  expected = (-1).^(cs*mm);
  printf('  CS_phase=%d  residual %.1e  max off-diagonal %.1e  diagonal is %s\n', ...
     cs, max(max(abs(Bm*T - Ys))), offdiag, ...
     merge(max(abs(diag(T)'-expected)) < 1e-5, ...
           merge(cs==0,'the identity','(-1)^m'), 'SOMETHING ELSE'));
end
" 2>/dev/null | grep CS_phase

echo "=== 2. where MRtrix's sh2peaks finds a fODF written in SMI's convention"
echo "    true axis:  0.3010 -0.5017  0.8127"
for cs in 0 1; do
  sh2peaks -quiet -num 1 $W/fod_cs$cs.mih $W/pk$cs.mih -force
  v=$(mrconvert -quiet $W/pk$cs.mih -coord 0 0 -coord 1 0 -coord 2 0 - | mrdump -quiet - | tr '\n' ' ')
  echo "    CS_phase=$cs peak (unnormalised): $v"
done
python3 - <<'PY'
import subprocess, numpy as np
n = np.array([0.30, -0.50, 0.81]); n /= np.linalg.norm(n)
for cs in (0, 1):
    out = subprocess.run(
        'mrconvert -quiet mrtrix/basis_check/pk%d.mih -coord 0 0 -coord 1 0 '
        '-coord 2 0 - | mrdump -quiet -' % cs, shell=True, capture_output=True)
    v = np.array([float(x) for x in out.stdout.split()])
    v = v / np.linalg.norm(v)
    ang = np.degrees(np.arccos(np.clip(abs(v @ n), 0, 1)))
    print(f'    CS_phase={cs}: {np.round(v,4)}  ->  {ang:6.2f} deg from the truth')
PY

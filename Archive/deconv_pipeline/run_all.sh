#!/bin/bash
# Everything in Reports/REPORT_SMI_deconvolution_MonteCarlo.md, from an empty data/.
#
# MATLAB/Octave does SMI and writes the images; MRtrix3 does CSD, MSMT-CSD and
# every peak extraction; Python does the bookkeeping and the figures.
#
#   ./run_all.sh            full run, 10,000 realisations per condition
#   NREP=200 ./run_all.sh   a quick version that exercises every step
set -e
cd "$(dirname "$0")"
NREP=${NREP:-10000}
NPROC=${NPROC:-4}
OCT="${OCTAVE:-octave-cli} --no-gui -q --eval"
P="$(pwd)/oct_path.m"
export OMP_NUM_THREADS=1

echo "== protocol"
python3 setup_protocol.py

echo "== phantom and response functions"
$OCT "run('$P'); gen_phantom(30,'p30');"
$OCT "run('$P'); dump_bases;"
./run_mrtrix.sh responses

echo "== SH basis check against MRtrix"
./check_mrtrix_basis.sh

echo "== Monte Carlo arms"
for s in 50 30 20 10 5; do
  $OCT "run('$P'); gen_montecarlo($s,$NREP,'snr$s');" &
  while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 10; done
done
$OCT "run('$P'); gen_montecarlo(1e4,$((NREP/20+50)),'nf');"
wait

for t in snr50 snr30 snr20 snr10 snr5 nf; do
  ./run_mrtrix.sh fit $t
done

echo "== tables and figures"
python3 tables.py nf snr50:50 snr30:30 snr20:20 snr10:10 snr5:5 \
    | tee ../Reports/deconv_tables.md
python3 figure_mc.py ../Figures/fodf_deconv_montecarlo.png \
    snr50:50 snr30:30 snr20:20 snr10:10 snr5:5
python3 figure_response.py p30 ../Figures/fodf_response_shview.png

echo "== the non-negativity sweep (report section 6)"
$OCT "run('$P'); sweep_nonneg(30,$((NREP/5)),'sw30');"
./run_mrtrix.sh sweep sw30 6
python3 score_sweep.py sw30 30

echo "== the exact-response control (report section 7)"
./run_mrtrix.sh control snr50
python3 -c "
import score_mrtrix as sm
sm.ARMS = [('SSST-CSD exact','csdXfod'), ('MSMT-CSD exact','msmtXfod')]
sm.report('snr50', 50)"

#!/bin/bash
# Everything in REPORT_SMI_deconvolution_MonteCarlo.md, from an empty data/.
#
# The five SNR arms are independent; NPROC controls how many run at once and
# how many workers the MSMT-CSD fits use. About two hours on four cores.
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

echo "== protocol and evaluation sphere"
python3 setup_protocol.py

echo "== phantom and response functions"
$OCT "run('$P'); gen_phantom(30,'p30');"
$OCT "run('$P'); dump_bases;"
python3 dhollander.py p30

echo "== Monte Carlo arms"
for s in 50 30 20 10 5; do
  $OCT "run('$P'); gen_montecarlo($s,$NREP,'snr$s');" &
  while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 10; done
done
$OCT "run('$P'); gen_montecarlo(1e4,$((NREP/20+50)),'nf');"
wait

for t in snr50 snr30 snr20 snr10 snr5 nf; do
  python3 run_csd_par.py $t p30 $NPROC
done

echo "== checks, tables and figures"
python3 check_conventions.py snr30
python3 tables.py nf snr50:50 snr30:30 snr20:20 snr10:10 snr5:5 \
    | tee ../deconv_tables.md
python3 figure_mc.py ../fodf_deconv_montecarlo.png \
    snr50:50 snr30:30 snr20:20 snr10:10 snr5:5
python3 figure_response.py p30 ../fodf_response_shview.png

echo "== the non-negativity sweep and the exact-response control"
$OCT "run('$P'); sweep_nonneg(30,$((NREP/10)),'sw30');"
python3 score_sweep.py sw30 30
python3 control_exact.py snr30 500
python3 score_control.py snr30X 30

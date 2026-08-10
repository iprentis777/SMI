#!/bin/bash
# CSD and MSMT-CSD, by MRtrix3 itself.
#
# Nothing in this file is a reimplementation. `dwi2response`, `dwi2fod` and
# `sh2peaks` are the MRtrix3 binaries; the only thing the rest of the package
# does is write the images they read (mrtrix_io.m) and read the images they
# write. The SMI fODF goes through the SAME `sh2peaks`, so the peak finder is
# MRtrix's for every arm including SMI's.
#
#   ./run_mrtrix.sh responses          estimate the response functions
#   ./run_mrtrix.sh fit <tag>          deconvolve one Monte Carlo arm
#
set -e
cd "$(dirname "$0")"
MD=mrtrix
NUMPEAKS=4
LMAX=6

case "$1" in

responses)
  # dwi2response dhollander: unsupervised 3-tissue estimation, the MRtrix
  # default for MSMT-CSD. -erode 0 because the phantom has no brain edge for
  # the erosion to remove.
  rm -rf $MD/scratch_dh
  dwi2response dhollander $MD/phantom_p30.mih \
      $MD/resp_wm.txt $MD/resp_gm.txt $MD/resp_csf.txt \
      -mask $MD/phantom_mask.mih -erode 0 -lmax 0,$LMAX,$LMAX,$LMAX \
      -voxels $MD/dh_voxels.mih -scratch $MD/scratch_dh -force
  # dwi2response tournier: the single-shell pairing for dwi2fod csd, on the
  # highest shell, which is where single-shell CSD is normally run.
  # dwi2response is a Python script and cannot read from a pipe, so the
  # single-shell subset goes through a file.
  dwiextract $MD/phantom_p30.mih -shells 0,3000 $MD/phantom_b3.mih -force -quiet
  rm -rf $MD/scratch_t $MD/scratch_fa
  dwi2response tournier $MD/phantom_b3.mih $MD/resp_wm_tournier.txt \
      -mask $MD/phantom_mask.mih -lmax $LMAX \
      -voxels $MD/tournier_voxels.mih -scratch $MD/scratch_t -force
  # dwi2response fa, the old selector, for comparison only
  dwi2response fa $MD/phantom_b3.mih $MD/resp_wm_fa.txt \
      -mask $MD/phantom_mask.mih -lmax $LMAX -scratch $MD/scratch_fa -force
  echo "--- WM responses (rows = shells, columns = even l)"
  echo "dhollander:"; cat $MD/resp_wm.txt
  echo "tournier:";   cat $MD/resp_wm_tournier.txt
  echo "fa:";         cat $MD/resp_wm_fa.txt
  ;;

fit)
  T=$2
  dwi2fod msmt_csd $MD/mc_$T.mih \
      $MD/resp_wm.txt  $MD/msmtfod_$T.mih \
      $MD/resp_gm.txt  $MD/msmtgm_$T.mih \
      $MD/resp_csf.txt $MD/msmtcsf_$T.mih \
      -mask $MD/mask_$T.mih -lmax $LMAX,0,0 -force -quiet
  dwiextract $MD/mc_$T.mih -shells 0,3000 $MD/mc_b3_$T.mih -force -quiet
  dwi2fod csd $MD/mc_b3_$T.mih $MD/resp_wm_tournier.txt $MD/csdfod_$T.mih \
      -mask $MD/mask_$T.mih -lmax $LMAX -force -quiet
  for f in smifod csdfod msmtfod gtfod; do
    sh2peaks $MD/${f}_$T.mih $MD/pk_${f}_$T.mih \
        -num $NUMPEAKS -mask $MD/mask_$T.mih -force -quiet
  done
  echo "fit $T done"
  ;;

control)
  # the same two methods, given the exact response of the kernel that generated
  # the data instead of an estimated one
  T=$2
  python3 write_exact_response.py $MD
  dwi2fod msmt_csd $MD/mc_$T.mih \
      $MD/resp_wm_exact.txt  $MD/msmtXfod_$T.mih \
      $MD/resp_gm_exact.txt  $MD/msmtXgm_$T.mih \
      $MD/resp_csf_exact.txt $MD/msmtXcsf_$T.mih \
      -mask $MD/mask_$T.mih -lmax $LMAX,0,0 -force -quiet
  dwi2fod csd $MD/mc_b3_$T.mih $MD/resp_wm_exact_b3.txt $MD/csdXfod_$T.mih \
      -mask $MD/mask_$T.mih -lmax $LMAX -force -quiet
  for f in csdXfod msmtXfod; do
    sh2peaks $MD/${f}_$T.mih $MD/pk_${f}_$T.mih \
        -num $NUMPEAKS -mask $MD/mask_$T.mih -force -quiet
  done
  echo "control $T done"
  ;;

sweep)
  # the lambda_nonneg sweep: same sh2peaks, so its numbers sit on the same
  # scale as the four-arm comparison
  T=$2
  N=$3
  for i in $(seq 1 $N); do
    sh2peaks $MD/sweep${i}fod_$T.mih $MD/pk_sweep${i}fod_$T.mih \
        -num $NUMPEAKS -mask $MD/mask_$T.mih -force -quiet
  done
  sh2peaks $MD/gtfod_$T.mih $MD/pk_gtfod_$T.mih \
      -num $NUMPEAKS -mask $MD/mask_$T.mih -force -quiet
  echo "sweep peaks $T done"
  ;;

*)
  echo "usage: $0 responses | fit <tag> | control <tag> | sweep <tag> <n>"; exit 1;;
esac

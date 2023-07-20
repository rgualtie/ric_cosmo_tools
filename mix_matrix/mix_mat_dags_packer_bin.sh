#!/bin/bash

# make the files with
# /big_scratch/rgualtie/condor_logs/$s\_0$b\_s10$i\_idx0$i/$s\_0$b\_s10$i\_idx0$i\_$fr\GHz_clean.dag
# Arguments are frequency, bin (3 digits) and spectral type: TT EE BB

fr=$1
b=$2
s=$3

for i in {0..99};do 
    python src/mix_matrix_parallel_2.1.py --nside 2048 --infile /sptgrid/user/rgualtie/mix_matrix/input_spectra/spectra_$s\_0$b.npy --band $fr --simidx $i
    echo "$i: Waiting 1sec.."
    sleep 1
done

files=/big_scratch/rgualtie/condor_logs/$s\_0$b\_s1*_idx*/$s\_0$b\_s1*_idx*_$fr\GHz_clean.dag
echo $files
for f in $files;do 
    (cat "${f}"; echo) >> /big_scratch/rgualtie/condor_logs/mix_mat/mix_mat_$fr\GHz_$s\_$b.dag
done
condor_submit_dag -update_submit /big_scratch/rgualtie/condor_logs/mix_mat/mix_mat_$fr\GHz_$s\_$b.dag

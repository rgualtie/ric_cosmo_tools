#!/bin/bash

# make the files with
# for fr in 90 150 220;do for s in TT EE BB;do for b in {008..508..25};do for i in {0..99};do python src/mix_matrix_parallel_2.1.py --nside 512 --infile /sptgrid/user/rgualtie/mix_matrix/input_spectra/spectra_$s\_0$b.npy --band $fr --simidx $i;done;done;done;done
# the dags will look like:
# /big_scratch/rgualtie/condor_logs/$s\_0$b\_s10$i\_idx0$i/$s\_0$b\_s10$i\_idx0$i\_$fr\GHz_clean.dag

for s in 90;do for i in {0..99};do python src/mix_matrix_parallel_2.1.py --nside 512 --infile /sptgrid/user/rgualtie/mix_matrix/input_spectra/spectra_TT_0483.npy --band $s --simidx $i;done;done

# the index i represents the decades 
for i in {00..09};do
#for i in 00;do
# Fr for frequencies
#    for fr in 90 150 220;do
    for fr in 90;do
# s for spectral type
#        for s in TT EE BB;do
        for s in TT;do
# b for bins
#            for b in {008.508..25};do
            for b in 483;do
                files=/big_scratch/rgualtie/condor_logs/$s\_0$b\_s1$i*_idx$i*/$s\_0$b\_s1$i*_idx$i*_$fr\GHz_clean.dag
                echo $files
                for f in $files;do 
                    (cat "${f}"; echo) >> /big_scratch/rgualtie/condor_logs/mix_mat/mix_mat_$i\_$fr\GHz_$s\_$b.dag
                done
                condor_submit_dag -update_submit /big_scratch/rgualtie/condor_logs/mix_mat/mix_mat_$i\_$fr\GHz_$s\_$b.dag
                echo Wait 1min..
                sleep 60
            done
        done
    done
done

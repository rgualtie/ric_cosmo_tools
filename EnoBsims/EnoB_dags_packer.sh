#!/bin/bash
# Pass the frequency band as first argument
b=$1
start=$2
stop=$3

# Generate the dag files that will be packed in the next step
for ((i=$start; i<=$stop; i++));
    do python src/EnoB_parallel_2.1.py --band $b --simidx $i --cleanup
done
# Gather the files
files=/big_scratch/rgualtie/condor_logs/EnoBlens_s${start:0:1}*_idx*/EnoBlens_s${start:0:1}*_idx*_$b\GHz_clean.dag
echo $files
# Pack and submit
for f in $files;do 
    (cat "${f}"; echo) >> /big_scratch/rgualtie/condor_logs/EnoB/EnoB_$b\_$start:$stop.dag
done
condor_submit_dag -update_submit /big_scratch/rgualtie/condor_logs/EnoB/EnoB_$b\_$start:$stop.dag

#!/bin/bash
# Pass the frequency band as first argument
b=$1
# simulation range:
start=$2
stop=$3
binw=$4
# Generate the dag files that will be packed in the next step
for ((i=$start; i<=$stop; i++));do python src/EnoB_parallel_2.1.py --band $b --simidx $i --binw $binw --cleanup;done
# Gather the files
files=/big_scratch/$USER/condor_logs/EnoBlens_s1*_idx*/EnoBlens_s1*_idx*_$b\GHz_clean.dag
echo $files
# Pack and submit
for f in $files;do 
    (cat "${f}"; echo) >> /big_scratch/$USER/condor_logs/EnoB/EnoB_$b\GHz_$start:$stop.dag
done
condor_submit_dag -update_submit /big_scratch/$USER/condor_logs/EnoB/EnoB_$b\GHz_$start:$stop.dag

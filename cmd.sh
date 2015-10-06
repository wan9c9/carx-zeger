#!/bin/sh
carx_batch_variable_nObs=${1}
carx_batch_variable_cr_idx=${2}
export carx_batch_variable_nObs
export carx_batch_variable_cr_idx
echo $carx_batch_variable_nObs
echo $carx_batch_variable_cr_idx


mpirun -n 24 mpi-wrapper 1000 ". ./prepare.sh" ". ./simulation_single.sh " ". ./summarize.sh"  && \
mv summary.txt   summary_${carx_batch_variable_nObs}_${carx_batch_variable_cr_idx}.txt

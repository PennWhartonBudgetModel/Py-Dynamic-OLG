#!/bin/bash
#$ -q short.q
#$ -j y
#$ -o /dev/null


# Define deep parameter grid size and solution batch size
GRIDSIZE=40
BATCHSIZE=5

# Find number of batches
# (Should be less than the array job task limit of 75,000)
let NBATCHES=(${GRIDSIZE}**3+${BATCHSIZE}-1)/${BATCHSIZE}


# Submit deep parameter set generation job
qsub -N generate_ss_sets -q short.q -j y -o /dev/null -b y 'matlab -nojvm -nosplash -r "generate_ss_sets('${GRIDSIZE}','${BATCHSIZE}')"'


# Submit batch solution jobs, holding for deep parameter set generation job
qsub -N solve_ss_batch -hold_jid generate_ss_sets -t 1-${NBATCHES} -q short.q -j y -o /dev/null -b y 'matlab -nojvm -nosplash -r "solve_ss_batch(${SGE_TASK_ID})"'


# Submit target inversion job, holding for batch solution jobs
qsub -N invert_ss_targets -hold_jid solve_ss_batch -q short.q -j y -o /dev/null -b y 'matlab -nojvm -nosplash -r "generate_ss_inverter, invert_ss_targets"'
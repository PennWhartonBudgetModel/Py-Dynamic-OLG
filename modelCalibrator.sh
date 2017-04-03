#!/bin/bash
#$ -q short.q
#$ -j y
#$ -o /dev/null

# Build mex functions and define batches
matlab -nojvm -nosplash -r "mexBuilder.all(), modelCalibrator.define_batches()"

# Submit batch solution task array job
qsub -N solve_batch -t 1-$(ls -1 ./Batches | wc -l) -q short.q \
     -j y -o /dev/null \
     -b y 'matlab -nojvm -nosplash -r "modelCalibrator.solve_batch(${SGE_TASK_ID})"'

# Submit batch consolidation job, holding for batch solution task array job
qsub -N consolidate_batches -hold_jid solve_batch -q short.q \
     -j y -o /dev/null \
     -b y 'matlab -nojvm -nosplash -r "modelCalibrator.consolidate_batches(false)"'
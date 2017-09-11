#!/bin/bash
#$ -q short.q
#$ -j y
#$ -o /dev/null

# Build mex functions and define batches
matlab -nojvm -nosplash -r "mexBuilder.all(), modelCalibrator.define_batches()"

# Clear or create log directory
LOGDIR='./Logs'
rm -rf ${LOGDIR}
mkdir -p ${LOGDIR}

# Specify batch solution job execution queue, checking for presence of --aws flag
QUEUE=$([ $# -gt 0 ] && [ $1 = "--aws" ]               \
        && echo 'aws-ppi.q -l aws_ppi -P bepp_ppi_aws' \
        || echo 'short.q'                              )

# Submit batch solution task array job
qsub -N solve_batch -t 1-$(ls -1 ./Batches | wc -l) \
     -q ${QUEUE} \
     -j y -o ${LOGDIR}'/batch$TASK_ID.log'  \
     -b y 'matlab -nojvm -nosplash -r "modelCalibrator.solve_batch(${SGE_TASK_ID})"'

# Submit batch consolidation job, holding for batch solution task array job
qsub -N consolidate_batches -hold_jid solve_batch \
     -q short.q \
     -j y -o /dev/null \
     -b y 'matlab -nojvm -nosplash -r "modelCalibrator.consolidate_batches(false)"'

#!/bin/bash
#$ -N define
#$ -q short.q
#$ -j y
#$ -o /dev/null

# Build mex functions and define calibration points
matlab -nodesktop -nosplash -r "PathFinder.setToProductionMode(), MexBuilder.all(), ModelCalibrator.definePoints()"

# Clear or create log directory
LOGDIR='./Logs'
rm -rf ${LOGDIR}
mkdir -p ${LOGDIR}

# Specify calibration job execution queue, checking for presence of --aws flag ($1)
QUEUE=$([ $# -gt 0 ] && [ $1 = "--aws" ]               \
        && echo 'aws-ppi.q -l aws_ppi -P bepp_ppi_aws' \
        || echo 'short.q'                              )

# Submit calibration task array job
#   Note use of -nojvm to deactivate Java and hence deactivate Matlab parfor
qsub -N calibrate -t 1-$(ls -1 ./CalibrationPoints | wc -l) \
     -q ${QUEUE} \
     -j y -o ${LOGDIR}'/point$TASK_ID.log'  \
     -b y 'matlab -nojvm -nosplash -r "PathFinder.setToProductionMode(), ModelCalibrator.calibratePoint(${SGE_TASK_ID})"'

# Submit calibration point consolidation job, holding for calibration task array job
qsub -N consolidate -hold_jid calibrate \
     -q short.q \
     -j y -o ${LOGDIR}'/consolidate.log' \
     -b y 'matlab -nodesktop -nosplash -r "PathFinder.setToProductionMode(), ModelCalibrator.consolidatePoints()"'

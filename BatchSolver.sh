#!/bin/bash
#$ -N define
#$ -q short.q
#$ -j y
#$ -o /dev/null

# Build mex functions and define scenarios for specified batch ($1)
matlab -nodesktop -nosplash -r "PathFinder.setToProductionMode(), MexBuilder.all(), BatchSolver.defineScenarios($1)"

# Clear or create log directory
LOGDIR='./Logs'
rm -rf ${LOGDIR}
mkdir -p ${LOGDIR}

# Specify batch solution job execution queue, checking for presence of --aws flag ($2)
QUEUE=$([ $# -gt 1 ] && [ $2 = "--aws" ]               \
        && echo 'aws-ppi.q -l aws_ppi -P bepp_ppi_aws' \
        || echo 'short.q'                              )

# Submit batch solution task array job
# (Note use of -nojvm to deactivate Java and hence deactivate Matlab parfor)
qsub -N solve -t 1-$(ls -1 ./Scenarios | wc -l) \
     -q ${QUEUE} \
     -j y -o ${LOGDIR}'/scenario$TASK_ID.log' \
     -b y 'matlab -nojvm -nosplash -r "PathFinder.setToProductionMode(), BatchSolver.solve(${SGE_TASK_ID})"'

# Submit data series generation job, holding for batch solution task array job
qsub -N generate -hold_jid solve \
     -q short.q \
     -j y -o ${LOGDIR}'/generate.log' \
     -b y 'matlab -nodesktop -nosplash -r "PathFinder.setToProductionMode(), BatchSolver.checkTerminations(), BatchSolver.generateDataSeries($1)"'

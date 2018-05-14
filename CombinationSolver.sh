#!/bin/bash
#$ -N generate
#$ -q short.q
#$ -j y
#$ -o /dev/null

# Build mex functions and generate scenarios
matlab -nodesktop -nosplash -r "PathFinder.setToProductionMode(), MexBuilder.all(), CombinationSolver.generateScenarios();"

# Clear or create log directory
LOGDIR='./Logs'
rm -rf ${LOGDIR}
mkdir -p ${LOGDIR}

# Specify bundle solution job execution queue, checking for presence of --aws flag ($2)
QUEUE=$([ $# -gt 1 ] && [ $2 = "--aws" ]               \
        && echo 'aws-ppi.q -l aws_ppi -P bepp_ppi_aws' \
        || echo 'short.q'                              )

# Submit current policy solution task array job
#   Note use of -nojvm to deactivate Java and hence deactivate Matlab parfor
qsub -N currentpolicy \
     -t 1-$(ls -1 ./Scenarios/currentpolicy*.mat | wc -l) \
     -q ${QUEUE} \
     -j y -o ${LOGDIR}'/currentpolicy$TASK_ID.log' \
     -b y 'matlab -nojvm -nosplash -r "PathFinder.setToProductionMode(), CombinationSolver.solveCurrentPolicy(${SGE_TASK_ID})"'

# Submit counterfactual solution task array job, holding for current policy solution task array job
#   Note use of -nojvm to deactivate Java and hence deactivate Matlab parfor
qsub -N counterfactual -hold_jid currentpolicy \
     -t 1-$(ls -1 ./Scenarios/counterfactual*.mat | wc -l) \
     -q ${QUEUE} \
     -j y -o ${LOGDIR}'/counterfactual$TASK_ID.log' \
     -b y 'matlab -nojvm -nosplash -r "PathFinder.setToProductionMode(), CombinationSolver.solveCounterfactual(${SGE_TASK_ID})"'

# Submit series generation job, holding for counterfactual solution task array job
qsub -N series -hold_jid counterfactual \
     -q short.q \
     -j y -o ${LOGDIR}'/series.log' \
     -b y 'matlab -nodesktop -nosplash -r "PathFinder.setToProductionMode(), CombinationSolver.checkTerminations(), CombinationSolver.generateSeries()"'

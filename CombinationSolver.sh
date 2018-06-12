#!/bin/bash
#$ -N generate
#$ -q short.q
#$ -j y
#$ -o /dev/null

# Build mex functions and generate scenarios
matlab -nodesktop -nosplash -r "MexBuilder.all(), PathFinder.setToProductionMode(), CombinationSolver.generateScenarios();"

# Clear or create log directory
LOGDIR='./Logs'
rm -rf ${LOGDIR}
mkdir -p ${LOGDIR}

# Specify scenario job execution queue, checking for presence of --aws flag ($1)
QUEUE=$([ $# -gt 0 ] && [ $1 = "--aws" ]                               \
        && echo 'aws-ppi.q -l aws_ppi -P bepp_ppi_aws -pe openmp 6' \
        || echo 'short.q                              -pe openmp 2' )

# Submit current policy solution task array job
qsub -N currentpolicy \
     -t 1-$(ls -1 ./Scenarios/currentpolicy*.mat | wc -l) \
     -q ${QUEUE} \
     -j y -o ${LOGDIR}'/currentpolicy$TASK_ID.log' \
     -b y 'matlab -nodesktop -nosplash -r "parpool_hpcc(), PathFinder.setToProductionMode(), CombinationSolver.solveCurrentPolicy(${SGE_TASK_ID}), delete(gcp)"'

# Submit counterfactual solution task array job, holding for current policy solution task array job
qsub -N counterfactual -hold_jid currentpolicy \
     -t 1-$(ls -1 ./Scenarios/counterfactual*.mat | wc -l) \
     -q ${QUEUE} \
     -j y -o ${LOGDIR}'/counterfactual$TASK_ID.log' \
     -b y 'matlab -nodesktop -nosplash -r "parpool_hpcc(), PathFinder.setToProductionMode(), CombinationSolver.solveCounterfactual(${SGE_TASK_ID}), delete(gcp)"'

# Submit series generation job, holding for counterfactual solution task array job
qsub -N series -hold_jid counterfactual \
     -q short.q \
     -j y -o ${LOGDIR}'/series.log' \
     -b y 'matlab -nodesktop -nosplash -r "PathFinder.setToProductionMode(), CombinationSolver.checkTerminations(), CombinationSolver.generateSeries()"'

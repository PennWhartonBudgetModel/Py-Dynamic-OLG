#!/bin/bash
#$ -q short.q
#$ -j y
#$ -o /dev/null

# Build mex functions and define production runs
matlab -nojvm -nosplash -r "mexBuilder.all(), modelProducer.define_runs()"

# Clear or create log directory
LOGDIR='./Logs'
rm -rf ${LOGDIR}
mkdir -p ${LOGDIR}

# Specify production run job execution queue, checking for presence of --aws flag
QUEUE=$([ $# -gt 0 ] && [ $1 = "--aws" ]               \
        && echo 'aws-ppi.q -l aws_ppi -P bepp_ppi_aws' \
        || echo 'short.q'                              )

# Submit production run task array job
qsub -N run -t 1-$(ls -1 ./Runs | wc -l) \
     -q ${QUEUE} \
     -j y -o ${LOGDIR}'/run$TASK_ID.log' \
     -b y 'matlab -nojvm -nosplash -r "modelProducer.run(${SGE_TASK_ID})"'

# Submit results packaging job, holding for production run task array job
# (Note use of -nodesktop instead of -nojvm to activate Java for package_results method)
qsub -N package_results -hold_jid run \
     -q short.q \
     -j y -o ${LOGDIR}'/package_results.log' \
     -b y 'matlab -nodesktop -nosplash -r "modelProducer.check_terminations(), modelProducer.package_results()"'

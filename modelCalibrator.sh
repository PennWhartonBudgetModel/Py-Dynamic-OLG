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

# Submit inverter construction job, holding for batch solution task array job
# (Note use of -nodesktop instead of -nojvm to activate Java for plot_conditions method)
qsub -N construct_inverter -hold_jid solve_batch -q short.q \
     -j y -o /dev/null \
     -b y 'matlab -nodesktop -nosplash -r "modelCalibrator.construct_inverter(false), modelCalibrator.plot_conditions()"'
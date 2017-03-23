#!/bin/bash
#$ -q short.q
#$ -j y
#$ -o /dev/null


# Submit batch initialization job
qsub -N initialize_batches -q short.q -j y -o /dev/null -b y 'matlab -nojvm -nosplash -r "modelCalibrator.initialize_batches()"'


# Submit batch solution array job after batch initialization job completion
qsub -N submitter -hold_jid initialize_batches -q short.q -j y -o /dev/null -b y $'qsub -N solve_batch -t 1-$(ls -1 ./Batches | wc -l) -q short.q -j y -o /dev/null -b y \'matlab -nojvm -nosplash -r \"modelCalibrator.solve_batch(${SGE_TASK_ID})\"\''


# Submit inverter construction job, holding for batch solution array job completion
qsub -N construct_inverter -hold_jid submitter,solve_batch -q short.q -j y -o /dev/null -b y 'matlab -nojvm -nosplash -r "modelCalibrator.construct_inverter(false)"'
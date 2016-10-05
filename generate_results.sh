#!/bin/bash
#$ -q short.q
#$ -j y
#$ -o /dev/null


# Create empty log directory
LOGDIR='./Results/Logs'
rm -rf ${LOGDIR}
mkdir -p ${LOGDIR}


# Submit closed economy runs
for GCUT in 0.10 0.05 0.00 -0.05; do
  for PLAN in base trump clinton ryan; do
    qsub -N solve_closed -t 1-16 -q short.q -j y -o ${LOGDIR}'/closed_inddeep=$TASK_ID_plan='${PLAN}'_gcut='${GCUT}'.log' -b y 'matlab -nojvm -nosplash -r "solve_closed(inddeep_to_params(${SGE_TASK_ID}), '\'${PLAN}\'', '${GCUT}', false)"'
  done
done


# Submit open economy runs
# (Note that baseline runs are skipped since they are generated as required by the closed economy runs above)
for PLAN in trump clinton ryan; do
  qsub -N solve_open -t 1-16 -q short.q -j y -o ${LOGDIR}'/open_inddeep=$TASK_ID_plan='${PLAN}'.log' -b y 'matlab -nojvm -nosplash -r "solve_open(inddeep_to_params(${SGE_TASK_ID}), '\'${PLAN}\'', false)"'
done


# Process and package results
qsub -N package -hold_jid solve_open,solve_closed -q short.q -j y -o ${LOGDIR}'/package.log' -b y 'matlab -nojvm -nosplash -r "check_closed_convergence, generate_static_aggregates_all, package_results"'
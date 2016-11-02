#!/bin/bash
#$ -q short.q
#$ -j y
#$ -o /dev/null


# Create empty log directory
LOGDIR='./Logs'
rm -rf ${LOGDIR}
mkdir -p ${LOGDIR}


# Submit closed economy runs
# (Steady state and open economy runs executed as dependencies)
for GCUT in +0.00 +0.10 +0.05 -0.05; do
  for PLAN in base trump clinton ryan; do
    
    # Skip non-zero gcut baseline runs
    if [[ (${PLAN} == base) && (${GCUT} != +0.00) ]]; then continue; fi
    
    qsub -N solve_closed -t 1-16 -q short.q \
         -j y -o ${LOGDIR}'/closed_inddeep=$TASK_ID_plan='${PLAN}'_gcut='${GCUT}'.log' \
         -b y 'matlab -nojvm -nosplash -r "solve_closed(inddeep_to_params(${SGE_TASK_ID}), '\'${PLAN}\'', '${GCUT}', false)"'
		
  done
done


# Process and package results
qsub -N package -hold_jid solve_closed -q short.q \
     -j y -o ${LOGDIR}'/package.log' \
     -b y 'matlab -nojvm -nosplash -r "check_closed_convergence, generate_static_aggregates_closed, package_results"'
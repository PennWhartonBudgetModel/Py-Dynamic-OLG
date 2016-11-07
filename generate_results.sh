#!/bin/bash
#$ -q short.q
#$ -j y
#$ -o /dev/null


# Build mex files
./build_all.sh


# Create empty log directory
LOGDIR='./Logs'
rm -rf ${LOGDIR}
mkdir -p ${LOGDIR}


# Submit closed economy baseline runs
# (Steady state and open economy baseline runs performed as dependencies)
qsub -N baseline -t 1-16 -q short.q \
     -j y -o ${LOGDIR}'/closed_inddeep=$TASK_ID_base.log' \
     -b y 'sleep $[(${SGE_TASK_ID}-1)*5]; matlab -nodesktop -nosplash -r "pool = parpool(12), solve_closed(inddeep_to_params(${SGE_TASK_ID}), '\''base'\'', +0.00, false), delete(pool)"'


# Submit closed economy counterfactual runs
# (Corresponding open economy counterfactual runs performed as dependencies)
for GCUT in +0.10 +0.05 +0.00 -0.05; do
  for PLAN in trump clinton ryan; do
    
		# (Note that parallelization across cohorts is disabled with -nojvm)
    qsub -N counterfactual -hold_jid baseline -t 1-16 -q short.q \
         -j y -o ${LOGDIR}'/closed_inddeep=$TASK_ID_plan='${PLAN}'_gcut='${GCUT}'.log' \
         -b y 'matlab -nojvm -nosplash -r "solve_closed(inddeep_to_params(${SGE_TASK_ID}), '\'${PLAN}\'', '${GCUT}', false)"'
		
  done
done


# Process and package results
qsub -N package -hold_jid counterfactual -q short.q \
     -j y -o ${LOGDIR}'/package.log' \
     -b y 'matlab -nodesktop -nosplash -r "pool = parpool(12), check_closed_convergence, generate_static_aggregates_closed, package_results, delete(pool)"'
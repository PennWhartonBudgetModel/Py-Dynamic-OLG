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


# Submit closed economy counterfactual runs
# (Steady state, open economy, and baseline runs all performed as dependencies)
for GCUT in +0.10 +0.05 +0.00 -0.05; do
  for PLAN in trump clinton ryan; do
    
    qsub -N counterfactual -t 1-16 -q short.q \
         -j y -o ${LOGDIR}'/closed_inddeep=$TASK_ID_plan='${PLAN}'_gcut='${GCUT}'.log' \
         -b y 'matlab -nojvm -nosplash -r "dynamicSolver.closed( get_basedef(${SGE_TASK_ID}), struct('\''plan'\'', '\'${PLAN}\'', '\''gcut'\'', '${GCUT}') )"'
		
  done
done


# Process and package results
# (Note use of -nodesktop instead of -nojvm to activate Java for package_results)
qsub -N package -hold_jid counterfactual -q short.q \
     -j y -o ${LOGDIR}'/package.log' \
     -b y 'matlab -nodesktop -nosplash -r "check_closed_convergence, package_results"'


#!/bin/bash
#$ -q short.q
#$ -j y
#$ -o /dev/null

# Copy all run log files
# Rem to use the correct LOGDIR
LOGDIR='./Logs'
cat run*.log > allruns.log

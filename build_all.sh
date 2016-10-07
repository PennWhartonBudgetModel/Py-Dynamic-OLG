#!/bin/bash
#$ -q short.q
#$ -j y
#$ -o /dev/null

matlab -nojvm -nosplash -r "mexBuilder.build_all"
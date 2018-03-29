#!/bin/bash
# Setup bash job headers

# load local environment

type1=full
type2=partial

for i in {1..10}
do
	sbatch single_R_run.sbatch $i $type1
	sbatch single_R_run.sbatch $i $type2
done
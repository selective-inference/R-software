#!/bin/bash
# Setup bash job headers

# load local environment

type1=knockoff
type2=knockoff+

for i in {1..10}
do
	sbatch single_R_run.sbatch $i $type1
	sbatch single_R_run.sbatch $i $type2
done
#!/bin/bash
# Setup bash job headers

# load local environment


for i in {1..10}
do
	sbatch single_R_run.sbatch $i
done
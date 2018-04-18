#!/bin/bash
# Setup bash job headers

# load local environment

rho1=0.2
rho2=0.4
rho3=0.6

for i in {1..10}
do
	sbatch single_R_run.sbatch $i $rho1
	sbatch single_R_run.sbatch $i $rho2
	sbatch single_R_run.sbatch $i $rho3
done
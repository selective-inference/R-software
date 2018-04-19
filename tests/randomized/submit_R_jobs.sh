#!/bin/bash
# Setup bash job headers

# load local environment

type1=full
type2=partial
rho0=0
rho1=0.2
rho2=0.4
rho3=0.6

for i in {1..10}
do
	sbatch single_R_run.sbatch $i $type1 $rho0
	sbatch single_R_run.sbatch $i $type1 $rho1
	sbatch single_R_run.sbatch $i $type1 $rho2
	sbatch single_R_run.sbatch $i $type1 $rho3
	sbatch single_R_run.sbatch $i $type2 $rho0
	sbatch single_R_run.sbatch $i $type2 $rho1
	sbatch single_R_run.sbatch $i $type2 $rho2
	sbatch single_R_run.sbatch $i $type2 $rho3
done
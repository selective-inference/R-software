#!/bin/bash
# Setup bash job headers

# load local environment

method1=liu
method2=lee
method3=knockoff
method4=knockoff+
method5=randomized_full
method6=randomized_partial

for i in {1..1}
do
	sbatch single_R_run.sbatch $method1
	sbatch single_R_run.sbatch $method2
	sbatch single_R_run.sbatch $method3
	sbatch single_R_run.sbatch $method4
	sbatch single_R_run.sbatch $method5
	sbatch single_R_run.sbatch $method6
done
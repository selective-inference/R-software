setwd(getwd())
source("test_debiased_coverage.R")

args = commandArgs(trailingOnly=TRUE)
seed = round(as.numeric(args[1]))
rho = as.numeric(args[2])

outdir = paste("/scratch/users/jelenam/full/rho", toString(rho), "/", sep="")
outfile = file.path(outdir, paste(sep="","dl_result_", toString(seed), "_rho_", toString(rho), ".rds"))

test_debiased_coverage(seed=seed, outfile=outfile, rho=rho)


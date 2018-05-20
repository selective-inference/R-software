setwd(getwd())
source("test_knockoffs.R")

args = commandArgs(trailingOnly=TRUE)
seed = round(as.numeric(args[1]))
type=toString(args[2])
rho=as.numeric(args[3])

outdir = paste("/scratch/users/jelenam/full/rho", toString(rho), "/", sep="")
label=paste(type,"_result_", sep="")
outfile = file.path(outdir, paste(sep="",label, toString(seed), "_rho_", toString(rho), ".rds"))

test_knockoffs(seed=seed, outfile=outfile, method=type, rho=rho)


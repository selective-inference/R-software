setwd(getwd())
source("test_randomized.R")

args = commandArgs(trailingOnly=TRUE)
seed = round(as.numeric(args[1]))
type = toString(args[2])
rho = as.numeric(args[3])

outdir = paste("/scratch/users/jelenam/full/rho", toString(rho), "/", sep="")
label = paste("randomized_", type, "_result_", sep="")
outfile = file.path(outdir, paste(sep="",label, toString(seed), "_rho_", toString(rho), ".rds"))

test_randomized(seed=seed, outfile=outfile, type=type, rho=rho)


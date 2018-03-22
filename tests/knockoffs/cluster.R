setwd(getwd())
source("test_knockoffs.R")

args = commandArgs(trailingOnly=TRUE)
seed = round(as.numeric(args[1]))

outdir = "/scratch/users/jelenam/full/"
outfile = file.path(outdir, paste(sep="","knockoff_result_", toString(seed), ".rda"))

test_knockoffs(seed, outfile, method="knockoff")


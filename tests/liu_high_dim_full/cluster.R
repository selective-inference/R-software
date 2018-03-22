setwd(getwd())
source("test_liu_full.R")

args = commandArgs(trailingOnly=TRUE)
seed = round(as.numeric(args[1]))

outdir = "/scratch/users/jelenam/full/"
outfile = file.path(outdir, paste(sep="","liu_result_", toString(seed), ".rda"))

test_liu_full(seed, outfile)


setwd(getwd())
source("test_randomized_full.R")

args = commandArgs(trailingOnly=TRUE)
seed = round(as.numeric(args[1]))

outdir = "/scratch/users/jelenam/full/"
outfile = file.path(outdir, paste(sep="","randomized_result_", toString(seed), ".rda"))

test_randomized_full(seed, outfile)


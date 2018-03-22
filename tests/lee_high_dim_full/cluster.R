setwd(getwd())
source("test_lee_full.R")

args = commandArgs(trailingOnly=TRUE)
seed = round(as.numeric(args[1]))

outdir = "/scratch/users/jelenam/full/"
outfile = file.path(outdir, paste(sep="","lee_result_", toString(seed), ".rda"))

test_lee_full(seed, outfile)


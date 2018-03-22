setwd(getwd())
source("test_liu_full.R")

args = commandArgs(trailingOnly=TRUE)
seed = round(as.numeric(args[1]))

outdir = "/scratch/users/jelenam/full/"
outdir = paste(outdir, method, "_n", n, "p", p,"_", regen, sep='')
dir.create(outdir)
outfile = file.path(outdir, paste(sep="","list_result_", toString(seed), ".rda"))

test_liu_full(seed, outfile)


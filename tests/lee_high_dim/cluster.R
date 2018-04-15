setwd(getwd())
source("test_lee.R")

args = commandArgs(trailingOnly=TRUE)
seed = round(as.numeric(args[1]))
type = toString(args[2])

outdir = "/scratch/users/jelenam/full/"
label=paste("lee_", type,"_result_", sep="")
outfile = file.path(outdir, paste(sep="", label, toString(seed), ".rds"))

test_lee_full(seed=seed, outfile=outfile, type=type)


setwd(getwd())
source("test_randomized.R")

args = commandArgs(trailingOnly=TRUE)
seed = round(as.numeric(args[1]))
type = toString(args[2])
  
outdir = "/scratch/users/jelenam/full/"
label = paste("randomized_", type, "_result_", sep="")
outfile = file.path(outdir, paste(sep="",label, toString(seed), ".rds"))

test_randomized(seed, outfile, type=type)


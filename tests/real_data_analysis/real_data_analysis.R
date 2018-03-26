library(glmnet)
library(selectiveInference)
library(MASS)
library(knockoff)

setwd(getwd())
args = commandArgs(trailingOnly=TRUE)
method = toString(args[1])

outdir = "/scratch/users/jelenam/full/"
label = paste(method, "_results", sep="")
outfile = file.path(outdir, paste(sep="",label, ".rds"))

#setwd("/Users/Jelena/Dropbox/kevin/jelena/real_data")
set.seed(1)
loss="ls"
sigma_est = 0.62
lambda = 0.0352479
#lambda =  0.8*selectiveInference:::theoretical.lambda(X, loss, sigma_est) # "lambda" "0.0352479112219816"
#print(c("lambda", lambda))


liu_full = function(outfile){
  data=readRDS("real_data1.rds")
  X=data$X
  y=data$y
  n=nrow(X)
  penalty_factor = rep(1, ncol(X))
  
  soln = selectiveInference:::solve_problem_glmnet(X, y, lambda, penalty_factor=penalty_factor, loss=loss)
  PVS = selectiveInference:::inference_group_lasso(X, y, soln, groups=1:ncol(X), lambda=lambda, penalty_factor=penalty_factor, 
                                                 sigma_est, loss=loss, algo="glmnet", construct_ci = TRUE)
  saveRDS(list(active_vars=PVS$active_vars, 
               sel_intervals=PVS$sel_intervals, naive_intervals=PVS$naive_intervals,
               pvalues=PVS$pvalues, naive_pvalues=PVS$naive_pvalues), file=outfile)
  
  return(NULL)
}



lee_full = function(outfile){
  data=readRDS("real_data2.rds")
  X=data$X
  y=data$y
  n=nrow(X)

  lasso = glmnet(X, y, family=selectiveInference:::family_label(loss), alpha=1, standardize=FALSE, intercept=FALSE, thresh=1e-12)
  soln = as.numeric(coef(lasso,x=X,y=y, family=selectiveInference:::family_label(loss), s=lambda, exact=TRUE))[-1]
  PVS = selectiveInference:::fixedLassoInf(X,y,soln, intercept=FALSE, lambda*n, family=selectiveInference:::family_label(loss),
                                           type="full",sigma=sigma_est)
  
  abs_soln = abs(soln)
  beta_threshold = abs_soln[order(abs_soln,decreasing=TRUE)][length(PVS$pv)]
  active_vars = which(abs_soln>=beta_threshold)
  cat("nactive:", length(active_vars), "\n")
  cat("active vars:", active_vars, "\n")
  
  pvalues = PVS$pv
  sel_intervals = t(PVS$ci)
  
  saveRDS(list(active_vars=active_vars, sel_intervals=sel_intervals, pvalues=pvalues), file=outfile)
  
}


knockoff = function(method, outfile){
  if (method=="knockoff"){
    data=readRDS("real_data3.rds")
  } else if (method=="knockoff+"){
    data=readRDS("real_data4.rds")
  }
  
  X=data$X
  y=data$y
  
  offset=0
  if (method=="knockoff+"){
    offset=1
  }
  filter = knockoff.filter(X, y, fdr=q, offset=offset)
  saveRDS(list(active_vars = filter$rejected), file=outfile)
}


randomized = function(type, outfile){
  if (type=="full"){
    data=readRDS("real_data5.rds")
  } else if (type=="partial"){
    data=readRDS("real_data6.rds")
  }
  
  X=data$X
  y=data$y
  n=nrow(X)
  
  rand_lasso_soln = selectiveInference:::randomizedLasso(X, 
                                                         y, 
                                                         lambda*n, 
                                                         family=selectiveInference:::family_label(loss),
                                                         condition_subgrad=TRUE)
  
  full_targets=selectiveInference:::set.target(rand_lasso_soln, type=type, sigma_est=sigma_est)
  
  PVS = selectiveInference:::randomizedLassoInf(rand_lasso_soln,
                                                full_targets=full_targets,
                                                sampler = "norejection", #"adaptMCMC", #
                                                level=0.9, 
                                                burnin=1000, 
                                                nsample=10000)
  active_vars=rand_lasso_soln$active_set
  cat("active_vars:",active_vars,"\n")
  pvalues = PVS$pvalues
  sel_intervals = t(PVS$ci)  # matrix with two rows
  saveRDS(list(active_vars=active_vars, sel_intervals=sel_intervals, pvalues=pvalues), file=outfile)
  
  return(NULL)
}

if (method=="liu"){
  liu_full(outfile)
} else if (method=="lee"){
  lee_full(outfile)
} else if (method=="knockoff"){
  knockoff(method, outfile)
} else if (method=="knockoff+"){
  knockoff(method, outfile)
} else if (method=="randomized_full"){
  randomized("full", outfile)
} else if (method=="randomized_partial"){
  randomized("partial", outfile)
}






library(knockoff)
library(selectiveInference)


compute.knockoff = function(data, method, q=0.2, model.free=TRUE){
  
  X=data$X 
  y=data$y
  n = nrow(X)
  
  true.nonzero = which(data$beta!=0)
  s=length(true.nonzero)
  offset=0
  if (method=="knockoff_plus"){
    offset=1
  }
  
  if (model.free==TRUE){
    filter = knockoff.filter(X, y, fdr=q, offset=offset)
  } else{
    filter = knockoff.filter(X, y, knockoffs = create.fixed, fdr=q, offset=offset)
  }
  
  rejected = filter$selected
  nrejected = length(rejected)
  TP = length(intersect(rejected, true.nonzero))
  power = TP/s
  FDR = (nrejected-TP)/max(1, nrejected)
  
  return(list(rejected=rejected,
              power=power, 
              FDR=FDR, 
              nrejected = nrejected))
}



test_knockoffs = function(seed=1, outfile=NULL, method = "knockoff", loss="logit",
                          nrep=10, n=5000, p=100, s=20, rho=0.){
  
  snr = 5*sqrt(2*log(p)/n)
  #snr = 5*sqrt(2*log(p))
  
  set.seed(seed)

  FDR_sample=NULL
  power_sample=NULL
  
  for (i in 1:nrep){
    
    if (loss=="ls"){
      data = selectiveInference:::gaussian_instance(n=n, p=p, s=s, rho=rho, sigma=1, snr=snr)
    } else if (loss=="logit"){
      data = selectiveInference:::logistic_instance(n=n, p=p, s=s, rho=rho, snr=snr)
    }
    
    cat("true nonzero:", which(data$beta!=0), "\n")
    
    mc = compute.knockoff(data, method=method)
    FDR_sample=c(FDR_sample, mc$FDR)
    power_sample=c(power_sample, mc$power)
    
    if (length(FDR_sample)>0){
      print(c("FDR:", mean(FDR_sample)))
      print(c("power:", mean(power_sample)))
    }
  }
  
  if (is.null(outfile)){
    outfile=paste(method, ".rds", sep="")
  }
  
  saveRDS(list(FDR_sample=FDR_sample, power_sample=power_sample,
               n=n,p=p, s=s, snr=snr, rho=rho, method=method), file=outfile)
  
  return(list(FDR_sample=FDR_sample, power_sample=power_sample))
}

test_knockoffs()

larsInference=function(larsfit,x,y,sigma=NULL,nsteps=ncol(x)-1,alpha=0.1,trace=FALSE){
# 
# compute p-values and confidence intervals after each of the first nsteps of  the LAR algorithm.
# Inputs: larsfit- result of call to lars function, with type="lar"
#    x: matrix of predictors;    y: response
#    sigma= standard deviaion of noise; if null, it is estimated from the MSE of CV-estimated model
#  alpha= the desired (two-sided) confidence level
#   trace- flag to provide info along the way
# Outputs: pvalues.forw- exact (finite sample) pvalues for each predictor as it enters the model; based on
#   spacings approach of Taylor et al.
#    pvalues.remove- exact pvalues  for each predictor when it is removed from the model;
#    pvalues.forw.cov.unc- asymptotic  pvalues for each predictor as it enters the model,
#      based on the covariance test of Lockhart et al (2014), using Exp(1) null at each step
#     pvalues.forw.cov.cor- asymptotic  pvalues for each predictor as it enters the model,
#      based on the corrected covariance test of Taylor et al (2014),
#    selection.intervals- alpha-level confidence intervals for  each predictor at each LAR step
#      See Taylor et al (2014) for precise interpretation of these intervals
#     alpha- the desired (two-sided) confidence level used

n=length(y)
p=ncol(x)
  if (larsfit$type != "LAR") {
            stop("Call to Lars must use type='lar'")
        }

if(is.numeric(sigma)){sigmahat=sigma;
               if(trace){cat(c("Using fixed value sigma= ",sigma),fill=T)}
               }
if(is.null(sigma)){
   cvlars=cv.lars(x,y,type="lar",plot=F)
   khat=which.min(cvlars$cv)
   yhat=predict(larsfit,x,s=khat,mode="step",type="fit")$fit
   sigmahat=sqrt(sum((y-yhat)^2)/(n-khat-1))
if(trace){cat(c("Using estimated value sigma= ",sigmahat),fill=T)}
 }
ci=matrix(NA,ncol=2*nsteps,nrow=nsteps)
pv.remove=matrix(NA,ncol=nsteps,nrow=nsteps)
pv=pv.cov.cor=pv.cov.unc=rep(NA,nsteps)
nams=unlist(larsfit$act)[1:nsteps]
nams2=rep(NA,2*nsteps)
jj=0
for(j in 1:nsteps){
             jj=jj+1
    nams2[jj]=paste("  ",nams[j],"-lower",sep="")
              jj=jj+1
                   nams2[jj]=paste(nams[j],"-upper",sep="")
}
dimnames(ci)=list(paste("step ",1:nsteps),nams2)
dimnames(pv.remove)=list(paste("step ",1:nsteps),nams)
names(pv)=nams
names(pv.cov.unc)=nams
names(pv.cov.cor)=nams
 GL = genlasso(as.numeric(y),
                  x,
                  diag(1,ncol(x)),
                  maxsteps=min(nrow(x),ncol(x),nsteps+1),
                  approx=TRUE)

for(k in 1:nsteps){
if(trace){cat(c("step=",k),fill=T)}
larsinf.obj=lars_inference(x, y, GL, sigma=sigmahat, debugging=T,   ls_step=k, num_resid_contrasts=80)
 for(j in 1:k){
  a=selint(larsfit,j,k,alpha,sigma=sigmahat,trace=F,larsinf.obj=larsinf.obj)
  ci[k,1+(j-1)*2]=a$clow
  ci[k,2+(j-1)*2]=a$cup
  pv.remove[k,j]=a$pv
}
}
  pv=larsinf.obj$pv[,1]
 pv.cov.unc=larsinf.obj$pv[,2]
  pv.cov.cor=larsinf.obj$pv[,3]
weights=larsinf.obj$weights
names(pv)=nams
names(pv.cov.unc)=nams
names(pv)=nams
cat("",fill=T)
return(list(pvalues.forw=pv,pvalues.remove=pv.remove,pvalues.forw.cov.unc=pv.cov.unc, 
pvalues.forw.cov.cor=pv.cov.cor,selection.intervals=ci,alpha=alpha,sigmahat=sigmahat,weights=weights))
}


fixedCoxLassoInf=function(x, y, status,
                          beta, lambda,
                          alpha=.1, type=c("partial"),
                          tol.beta=1e-5, tol.kkt=0.1,
                          gridrange=c(-100,100), 
                          bits=NULL, verbose=FALSE,
                          this.call=NULL){

 checkargs.xy(x,y)
 if(is.null(status)) stop("Must supply `status' argument")
if( sum(status==0)+sum(status==1)!=length(y)) stop("status vector must have values 0 or 1")
  if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
  if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda") 
  checkargs.misc(beta=beta,lambda=lambda,alpha=alpha,
                 gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
  if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
    warning("Package Rmpfr is not installed, reverting to standard precision")
    bits = NULL
  }

    n=nrow(x)
    p=ncol(x)
     nvar=sum(beta!=0)
    pv=vlo=vup=sd=rep(NA, nvar)
     ci=tailarea=matrix(NA,nvar,2)


   m=beta!=0
vars=which(m)
if(sum(m)>0){
    bhat=beta[beta!=0] #penalized coefs just for active variables
    sign_bhat=sign(bhat)

 #check KKT
    
   aaa=coxph(Surv(y,status)~x[,m],init=bhat,iter.max=0) # this gives the Cox model at exactly bhat
                                                        # so when we compute gradient and score 
							# we are evaluating at the LASSO solution
							# naming of variables could be improved...
    res=residuals(aaa,type="score")
if(!is.matrix(res)) res=matrix(res,ncol=1)
scor=colSums(res)
    g=(scor+lambda*sign_bhat)/(2*lambda)
#    cat(c(g,lambda,tol.kkt),fill=T)
     if (any(abs(g) > 1+tol.kkt) )
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances)"))

# Hessian of partial likelihood at the LASSO solution    
MM=vcov(aaa)

bbar=(bhat+lambda*MM%*%sign_bhat)
A1=-(mydiag(sign_bhat))
b1= -(mydiag(sign_bhat)%*%MM)%*%sign_bhat*lambda

   temp=max(A1%*%bbar-b1)


# compute p-values

# JT: are we sure the signs of these are correctly handled?
# two sided p-values numerically agree with python but
# the one sided p-values are a bit off

    for(jj in 1:length(bbar)){
      vj=rep(0,length(bbar));vj[jj]=sign_bhat[jj]


      junk=TG.pvalue(bbar, A1, b1, vj,MM)

       pv[jj] = junk$pv
      vlo[jj]=junk$vlo
       vup[jj]=junk$vup
       sd[jj]=junk$sd

      junk2=TG.interval(bbar, A1, b1, vj, MM, alpha, flip=(sign_bhat[jj]==-1))
       ci[jj,]=junk2$int
       tailarea[jj,] = junk2$tailarea
     
  }
  # JT: these don't seem to be the real one-step estimators
    fit0=coxph(Surv(y,status)~x[,m])
      coef0=fit0$coef
      se0=sqrt(diag(fit0$var))
      zscore0=coef0/se0
    
   out = list(lambda=lambda,pv=pv,ci=ci,
    tailarea=tailarea,vlo=vlo,vup=vup,sd=sd,
    vars=vars,alpha=alpha,coef0=coef0,zscore0=zscore0,
    call=this.call)
  class(out) = "fixedCoxLassoInf" 
}
return(out)
}



print.fixedCoxLassoInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))
  
  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
  cat("",fill=T)
  tab = cbind(x$vars,
    round(x$coef0,3),
    round(x$zscore0,3),
    round(x$pv,3),round(x$ci,3))
  colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
  if (tailarea) {
    tab = cbind(tab,round(x$tailarea,3))
    colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
  }
  rownames(tab) = rep("",nrow(tab))
  print(tab)
 
  cat(sprintf("\nNote: coefficients shown are %s regression coefficients\n",
              ifelse(x$type=="partial","partial","full")))
  invisible()
}



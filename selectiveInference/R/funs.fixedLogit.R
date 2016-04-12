
fixedLogitLassoInf=function(x,y,beta,lambda,alpha=.1, type=c("partial"), tol.beta=1e-5, tol.kkt=0.1,
                     gridrange=c(-100,100), bits=NULL, verbose=FALSE,this.call=NULL){
    
    
  type = match.arg(type)
  checkargs.xy(x,y)
  if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
  if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda") 
  checkargs.misc(beta=beta,lambda=lambda,alpha=alpha,
                 gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
  if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
    warning("Package Rmpfr is not installed, reverting to standard precision")
    bits = NULL
  }

     
     n=length(y)
     p=ncol(x)
 # I assume that intcpt was used
     if(length(beta)!=p+1) stop("beta must be of length ncol(x)+1, that is, it should include an intercept")
     nvar=sum(beta[-1]!=0)
     pv=vlo=vup=sd=rep(NA, nvar)
     ci=tailarea=matrix(NA,nvar,2)
    
#do we need to worry about standardization?

#  obj = standardize(x,y,TRUE,FALSE)
  #  x = obj$x
  #  y = obj$y
  
 m=beta[-1]!=0  #active set
         
    bhat=c(beta[1],beta[-1][beta[-1]!=0]) # intcpt plus active vars
     s2=sign(bhat)
     lam2m=diag(c(0,rep(lambda,sum(m))))
 

    xxm=cbind(1,x[,m])
    
   etahat = xxm %*% bhat
   prhat = as.vector(exp(etahat) / (1 + exp(etahat)))
   ww=prhat*(1-prhat)
   w=diag(ww)
      
#check KKT
   z=etahat+(y-prhat)/ww
       g=  t(x)%*%w%*%(z-etahat)/lambda # negative gradient scaled by lambda
     if (any(abs(g) > 1+tol.kkt) )
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances)"))

 vars = which(abs(beta[-1]) > tol.beta / sqrt(colSums(x^2)))
  if(length(vars)==0){
      cat("Empty model",fill=T)
      return()
  }
  if (any(sign(g[vars]) != sign(beta[-1][vars])))
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances). You might try rerunning",
                  "glmnet with a lower setting of the",
                  "'thresh' parameter, for a more accurate convergence."))
  
 #constraints for active variables             
  MM=solve(t(xxm)%*%w%*%xxm)
  gm = c(0,-g[vars]*lambda) # gradient at LASSO solution, first entry is 0 because intercept is unpenalized
                            # at exact LASSO solution it should be s2[-1]
  dbeta = MM %*% gm

  # bbar=(bhat+lam2m%*%MM%*%s2)  # JT: this is wrong, shouldn't use sign of intercept anywhere...
  bbar = bhat - dbeta

  A1=-(mydiag(s2))[-1,]
  b1= (s2 * dbeta)[-1]

  tol.poly = 0.01 
  if (max((A1 %*% bbar)[-1] - b1) > tol.poly)
    stop(paste("Polyhedral constraints not satisfied; you must recompute beta",
               "more accurately. With glmnet, make sure to use exact=TRUE in coef(),",
               "and check whether the specified value of lambda is too small",
               "(beyond the grid of values visited by glmnet).",
               "You might also try rerunning glmnet with a lower setting of the",
               "'thresh' parameter, for a more accurate convergence."))


  
    for(jj in 1:sum(m)){
       vj=c(rep(0,sum(m)+1));vj[jj+1]=s2[jj+1]
      # compute p-values
      junk=mypoly.pval.lee(bbar,A1,b1,vj,MM)
      pv[jj] = junk$pv
 
   vlo[jj]=junk$vlo
   vup[jj]=junk$vup
       sd[jj]=junk$sd
  #  junk2=mypoly.int.lee(bbar[-1], A1, b1,vj,MM[-1,-1],alpha=.1)
     junk2=mypoly.int.lee(bbar,vj,vlo[jj],vup[jj],sd[jj],alpha=.1)

     ci[jj,]=junk2$int
     tailarea[jj,] = junk2$tailarea
   }

  # JT: these are not the one step estimators but they are close
  fit0=glm(y~x[,m],family="binomial")
  sfit0=summary(fit0)
      coef0=bbar[-1]        #fit0$coef[-1]
      se0=sqrt(diag(MM)[-1]) # sfit0$cov.scaled)[-1])
      zscore0=coef0/se0
      
  out = list(type=type,lambda=lambda,pv=pv,ci=ci,
    tailarea=tailarea,vlo=vlo,vup=vup,sd=sd,
    vars=vars,alpha=alpha,coef0=coef0,zscore0=zscore0,
    call=this.call,
    info.matrix=MM) # info.matrix is output just for debugging purposes at the moment
  class(out) = "fixedLogitLassoInf"
  return(out)

 }



print.fixedLogitLassoInf <- function(x, tailarea=TRUE, ...) {
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




# Lasso inference function (for fixed lambda). Note: here we are providing inference
# for the solution of
# min 1/2 || y - \beta_0 - X \beta ||_2^2 + \lambda || \beta ||_1

fixedLassoInf <- function(x, y, beta, lambda, family=c("gaussian","binomial","cox"),intercept=TRUE, add.targets=NULL, status=NULL,
sigma=NULL, alpha=0.1,
                     type=c("partial","full"), tol.beta=1e-5, tol.kkt=0.1,
                     gridrange=c(-100,100), bits=NULL, verbose=FALSE) {

  family = match.arg(family)
  this.call = match.call()
  type = match.arg(type)
  
  if(family=="binomial")  {
    if(type!="partial") stop("Only type= partial allowed with binomial family")
    out=fixedLogitLassoInf(x,y,beta,lambda,alpha=alpha, type="partial", tol.beta=tol.beta, tol.kkt=tol.kkt,
                           gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)
    return(out)
  }
  else if(family=="cox")  {
    if(type!="partial") stop("Only type= partial allowed with Cox family")
    out=fixedCoxLassoInf(x,y,status,beta,lambda,alpha=alpha, type="partial",tol.beta=tol.beta,
                         tol.kkt=tol.kkt, gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)
    return(out)
  }
  
  else{
    
    
    
    checkargs.xy(x,y)
    if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
    if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda")
    
    n = nrow(x)
    p = ncol(x)
    beta = as.numeric(beta)
    if (type == "full") {
      if (p > n) {
        # need intercept (if there is one) for debiased lasso
        hbeta = beta
        if (intercept == T) {
          if (length(beta) != p + 1) {
            stop("Since type='full', p > n, and intercept=TRUE, beta must have length equal to ncol(x)+1")
          }
          # remove intercept if included
          beta = beta[-1]
        } else if (length(beta) != p) {
          stop("Since family='gaussian', type='full' and intercept=FALSE, beta must have length equal to ncol(x)")
        }
      }
    } else if (length(beta) != p) {
      stop("Since family='gaussian' and type='partial', beta must have length equal to ncol(x)")
    }
    
    checkargs.misc(beta=beta,lambda=lambda,sigma=sigma,alpha=alpha,
                   gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
    if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
      warning("Package Rmpfr is not installed, reverting to standard precision")
      bits = NULL
    }
    
    if (!is.null(add.targets) && (!is.vector(add.targets)
                                  || !all(is.numeric(add.targets)) || !all(add.targets==floor(add.targets))
                                  || !all(add.targets >= 1 && add.targets <= p))) {
      stop("'add.targets' must be a vector of integers between 1 and p")
    }
    
    # If glmnet was run with an intercept term, center x and y
    if (intercept==TRUE) {
      obj = standardize(x,y,TRUE,FALSE)
      x = obj$x
      y = obj$y
    }
    
    # Check the KKT conditions
    g = t(x)%*%(y-x%*%beta) / lambda
    if (any(abs(g) > 1+tol.kkt * sqrt(sum(y^2))))
      warning(paste("Solution beta does not satisfy the KKT conditions",
                    "(to within specified tolerances)"))
    
    tol.coef = tol.beta * sqrt(n^2 / colSums(x^2))
    # print(tol.coef)
    vars = which(abs(beta) > tol.coef)
    # print(beta)
    # print(vars)
    if(length(vars)==0){
      cat("Empty model",fill=T)
      return()
    }
    if (any(sign(g[vars]) != sign(beta[vars])))
      warning(paste("Solution beta does not satisfy the KKT conditions",
                    "(to within specified tolerances). You might try rerunning",
                    "glmnet with a lower setting of the",
                    "'thresh' parameter, for a more accurate convergence."))
    
    # Get lasso polyhedral region, of form Gy >= u
    if (type == 'full' & p > n) out = fixedLasso.poly(x,y,beta,lambda,vars,inactive=TRUE)
    else out = fixedLasso.poly(x,y,beta,lambda,vars)
    G = out$G
    u = out$u
    
    # Check polyhedral region
    tol.poly = 0.01
    if (min(G %*% y - u) < -tol.poly * sqrt(sum(y^2)))
      stop(paste("Polyhedral constraints not satisfied; you must recompute beta",
                 "more accurately. With glmnet, make sure to use exact=TRUE in coef(),",
                 "and check whether the specified value of lambda is too small",
                 "(beyond the grid of values visited by glmnet).",
                 "You might also try rerunning glmnet with a lower setting of the",
                 "'thresh' parameter, for a more accurate convergence."))
    
    # Estimate sigma
    if (is.null(sigma)) {
      if (n >= 2*p) {
        oo = intercept
        sigma = sqrt(sum(lsfit(x,y,intercept=oo)$res^2)/(n-p-oo))
      }
      else {
        sigma = sd(y)
        warning(paste(sprintf("p > n/2, and sd(y) = %0.3f used as an estimate of sigma;",sigma),
                      "you may want to use the estimateSigma function"))
      }
    }
    
    # add additional targets for inference if provided
    if (!is.null(add.targets)) vars = sort(unique(c(vars,add.targets,recursive=T)))
    
    k = length(vars)
    pv = vlo = vup = numeric(k)
    vmat = matrix(0,k,n)
    ci = tailarea = matrix(0,k,2)
    sign = numeric(k)
      
    if (type=="full" & p > n) {
      if (intercept == T) {
        pp=p+1
        Xint <- cbind(rep(1,n),x)
        # indices of selected predictors
        S = c(1,vars + 1)
        notS = which(abs(beta) <= tol.coef) + 1
      } else {
        pp=p
        Xint <- x
        # indices of selected predictors
        S = vars
        notS = which(abs(beta) <= tol.coef)
      }
      
      
      XS = Xint[,S]
      hbetaS = hbeta[S]
      
      # Reorder so that active set S is first
      Xordered = Xint[,c(S,notS,recursive=T)]
      
      hsigma <- 1/n*(t(Xordered)%*%Xordered)
      hsigmaS <- 1/n*(t(XS)%*%XS) # hsigma[S,S]
      hsigmaSinv <- solve(hsigmaS) # pinv(hsigmaS)

      # Approximate inverse covariance matrix for when (n < p) from lasso_Inference.R
      htheta <- InverseLinfty(hsigma, n, length(S), verbose=FALSE)
      # htheta <- InverseLinfty(hsigma, n, verbose=FALSE)
      
      FS = rbind(diag(length(S)),matrix(0,pp-length(S),length(S)))
      GS = cbind(diag(length(S)),matrix(0,length(S),pp-length(S)))
      ithetasigma = (GS-(htheta%*%hsigma))
      # ithetasigma = (diag(pp) - (htheta%*%hsigma))
      
      M <- (((htheta%*%t(Xordered))+ithetasigma%*%FS%*%hsigmaSinv%*%t(XS))/n)
      # vector which is offset for testing debiased beta's
      null_value <- (((ithetasigma%*%FS%*%hsigmaSinv)%*%sign(hbetaS))*lambda/n)
      if (intercept == T) {
        M = M[-1,] # remove intercept row
        null_value = null_value[-1] # remove intercept element
      }
    } else if (type=="partial" || p > n) {
      xa = x[,vars,drop=F]
      M = pinv(crossprod(xa)) %*% t(xa)
      null_value = rep(0,k)
    } else {
      M = pinv(crossprod(x)) %*% t(x)
      M = M[vars,,drop=F]
      null_value = rep(0,k)
    }

  for (j in 1:k) {
    if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))

    vj = M[j,]
    mj = sqrt(sum(vj^2))
    vj = vj / mj        # Standardize (divide by norm of vj)
    sign[j] = sign(sum(vj*y))
    vj = sign[j] * vj

    limits.info = TG.limits(y, -G, -u, vj, Sigma=diag(rep(sigma^2, n)))
    a = TG.pvalue.base(limits.info, null_value=null_value[j], bits=bits)
    pv[j] = a$pv
    vlo[j] = a$vlo * mj # Unstandardize (mult by norm of vj)
    vup[j] = a$vup * mj # Unstandardize (mult by norm of vj)
    vmat[j,] = vj * mj * sign[j]  # Unstandardize (mult by norm of vj)

    a = TG.interval.base(limits.info, 
                         alpha=alpha,
                         gridrange=gridrange,
			 flip=(sign[j]==-1),
                         bits=bits)
    ci[j,] = (a$int-null_value[j]) * mj # Unstandardize (mult by norm of vj)
    tailarea[j,] = a$tailarea
  }

  out = list(type=type,lambda=lambda,pv=pv,ci=ci,
    tailarea=tailarea,vlo=vlo,vup=vup,vmat=vmat,y=y,
    vars=vars,sign=sign,sigma=sigma,alpha=alpha,
    sd=sigma*sqrt(rowSums(vmat^2)),
    coef0=vmat%*%y,
    call=this.call)
  class(out) = "fixedLassoInf"
  return(out)
}
}

#############################


fixedLasso.poly=
  function(x, y, beta, lambda, a, inactive = FALSE) {
    xa = x[,a,drop=F]
    xac = x[,!a,drop=F]
    xai = pinv(crossprod(xa))
    xap = xai %*% t(xa)
    za = sign(beta[a])
    if (length(za)>1) dz = diag(za)
    if (length(za)==1) dz = matrix(za,1,1)
    
    if (inactive) {
      P = diag(1,nrow(xa)) - xa %*% xap
      
      G = -rbind(
        1/lambda * t(xac) %*% P,
        -1/lambda * t(xac) %*% P,
        -dz %*% xap
      )
      lambda2=lambda
      if(length(lambda)>1) lambda2=lambda[a]
      u = -c(
        1 - t(xac) %*% t(xap) %*% za,
        1 + t(xac) %*% t(xap) %*% za,
        -lambda2 * dz %*% xai %*% za)
    } else {
      G = -rbind(
        #   1/lambda * t(xac) %*% P,
        # -1/lambda * t(xac) %*% P,
        -dz %*% xap
      )
      lambda2=lambda
      if(length(lambda)>1) lambda2=lambda[a]
      u = -c(
        #   1 - t(xac) %*% t(xap) %*% za,
        #   1 + t(xac) %*% t(xap) %*% za,
        -lambda2 * dz %*% xai %*% za)
    }
    
    return(list(G=G,u=u))
  }

##############################

### Functions borrowed and slightly modified from lasso_inference.R

## Approximates inverse covariance matrix theta
InverseLinfty <- function(sigma, n, e, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
  # InverseLinfty <- function(sigma, n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
    isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, e, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:e) {
    if ((i %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }
        }
      }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M)
}

InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=50, threshold=1e-2 ) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);
  
  if (mu >= mu0){
    beta[i] <- (1-mu0)/sigma[i,i];
    returnlist <- list("optsol" = beta, "iter" = 0);
    return(returnlist);
  }
  
  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta[i] <- (1-mu0)/sigma[i,i];
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta;
  
  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){
    
    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i)
        v <- v+1;
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      if (oldval != beta[j]){
        vs <- vs + (oldval-beta[j])*sigma.tilde[,j];
      }
    }
    
    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      if (iter>10)
        vs <- -sigma.tilde%*%beta;
    }
  }
  
  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}

SoftThreshold <- function( x, lambda ) {
  #
  # Standard soft thresholding
  #
  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}

##############################

print.fixedLassoInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))

  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
  cat("",fill=T)
  tab = cbind(x$vars,
    round(x$coef0,3),
    round(x$coef0 / x$sd,3),
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

#estimateLambda <- function(x, sigma, nsamp=1000){
#  checkargs.xy(x,rep(0,nrow(x)))
#  if(nsamp < 10) stop("More Monte Carlo samples required for estimation")
#  if (length(sigma)!=1) stop("sigma should be a number > 0")
 # if (sigma<=0) stop("sigma should be a number > 0")

 # n = nrow(x)
 # eps = sigma*matrix(rnorm(nsamp*n),n,nsamp)
 # lambda = 2*mean(apply(t(x)%*%eps,2,max))
 # return(lambda)
#}


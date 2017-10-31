# First part is only functions from the old code. At the bottom is
# the bit of code that actually compares the old vs new code

######################################################

### Old code (using Adel's R code)

## Approximates inverse covariance matrix theta
InverseLinfty <- function(sigma, n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-10, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:p) {
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

InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=50, threshold=1e-10) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);
  
  #if (mu >= mu0){
  #  beta[i] <- (1-mu0)/sigma[i,i];
  #  returnlist <- list("optsol" = beta, "iter" = 0);
  #  return(returnlist);
  #}
  
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
      #if (iter>10)
      #  vs <- -sigma.tilde%*%beta;
    }

    # print(c(iter, maxiter, diff.norm2, threshold * last.norm2, threshold, mu))

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


### Functions borrowed from selective Inference (only fixedLassoInf and fixedLasso.poly are modified)

# Special linear time order function, works only when x
# is a scrambled vector of integers.

Order <- function(x) {
  n = length(x)
  o = numeric(n)
  o[x] = Seq(1,n)
  return(o)
}

# Returns a sequence of integers from a to b if a <= b,
# otherwise nothing. You have no idea how important this
# function is...

Seq <- function(a, b, ...) {
  if (a<=b) return(seq(a,b,...))
  else return(numeric(0))
}

# Returns the sign of x, with Sign(0)=1.

Sign <- function(x) {
  return(-1+2*(x>=0))
}

##############################

# Centering and scaling convenience function

standardize <- function(x, y, intercept, normalize) {
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  if (intercept) {
    bx = colMeans(x)
    by = mean(y)
    x = scale(x,bx,FALSE)
    y = y-mean(y)
  } else {
    bx = rep(0,p)
    by = 0
  }
  if (normalize) {
    sx = sqrt(colSums(x^2))
    x = scale(x,FALSE,sx)
  } else {
    sx = rep(1,p)
  }
  
  return(list(x=x,y=y,bx=bx,by=by,sx=sx))
}

##############################

# Interpolation function to get coefficients

coef.interpolate <- function(betas, s, knots, dec=TRUE) {
  # Sort the s values
  o = order(s,dec=dec)
  s = s[o]
  
  k = length(s)
  mat = matrix(rep(knots,each=k),nrow=k)
  if (dec) b = s >= mat
  else b = s <= mat
  blo = max.col(b,ties.method="first")
  bhi = pmax(blo-1,1)
  
  i = bhi==blo
  p = numeric(k)
  p[i] = 0
  p[!i] = ((s-knots[blo])/(knots[bhi]-knots[blo]))[!i]
  
  beta = t((1-p)*t(betas[,blo,drop=FALSE]) + p*t(betas[,bhi,drop=FALSE]))
  colnames(beta) = as.character(round(s,3))
  rownames(beta) = NULL
  
  # Return in original order
  o = order(o)
  return(beta[,o,drop=FALSE])
}

##############################

checkargs.xy <- function(x, y) {
  if (missing(x)) stop("x is missing")
  if (is.null(x) || !is.matrix(x)) stop("x must be a matrix")
  if (missing(y)) stop("y is missing")
  if (is.null(y) || !is.numeric(y)) stop("y must be numeric")
  if (ncol(x) == 0) stop("There must be at least one predictor [must have ncol(x) > 0]")
  if (checkcols(x)) stop("x cannot have duplicate columns")
  if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 0]")
  if (length(y)!=nrow(x)) stop("Dimensions don't match [length(y) != nrow(x)]")
}

checkargs.misc <- function(sigma=NULL, alpha=NULL, k=NULL,
                           gridrange=NULL, gridpts=NULL, griddepth=NULL,
                           mult=NULL, ntimes=NULL,
                           beta=NULL, lambda=NULL, tol.beta=NULL, tol.kkt=NULL,
                           bh.q=NULL) {
  
  if (!is.null(sigma) && sigma <= 0) stop("sigma must be > 0")
  if (!is.null(lambda) && lambda < 0) stop("lambda must be >= 0")
  if (!is.null(alpha) && (alpha <= 0 || alpha >= 1)) stop("alpha must be between 0 and 1")
  if (!is.null(k) && length(k) != 1) stop("k must be a single number")
  if (!is.null(k) && (k < 1 || k != floor(k))) stop("k must be an integer >= 1")
  if (!is.null(gridrange) && (length(gridrange) != 2 || gridrange[1] > gridrange[2]))
    stop("gridrange must be an interval of the form c(a,b) with a <= b")
  if (!is.null(gridpts) && (gridpts < 20 || gridpts != round(gridpts)))
    stop("gridpts must be an integer >= 20")
  if (!is.null(griddepth) && (griddepth > 10 || griddepth != round(griddepth)))
    stop("griddepth must be an integer <= 10")
  if (!is.null(mult) && mult < 0) stop("mult must be >= 0")
  if (!is.null(ntimes) && (ntimes <= 0 || ntimes != round(ntimes)))
    stop("ntimes must be an integer > 0")
  if (!is.null(beta) && sum(beta!=0)==0) stop("Value of lambda too large, beta is zero")
  # if (!is.null(lambda) && length(lambda) != 1) stop("lambda must be a single number")
  if (!is.null(lambda) && length(lambda) != 1 && length(lambda) != length(beta)) stop("lambda must be a single number or equal to the length of beta")
  if (!is.null(lambda) && lambda < 0) stop("lambda must be >=0")
  if (!is.null(tol.beta) && tol.beta <= 0) stop("tol.beta must be > 0")
  if (!is.null(tol.kkt) && tol.kkt <= 0) stop("tol.kkt must be > 0")
}

# Make sure that no two columms of A are the same
# (this works with probability one).

checkcols <- function(A) {
  b = rnorm(nrow(A))
  a = sort(t(A)%*%b)
  if (any(diff(a)==0)) return(TRUE)
  return(FALSE)
}

estimateSigma <- function(x, y, intercept=TRUE, standardize=TRUE) {
  checkargs.xy(x,rep(0,nrow(x)))
  if(nrow(x)<10) stop("Number of observations must be at least 10 to run estimateSigma")
  cvfit=cv.glmnet(x,y,intercept=intercept,standardize=standardize)
  lamhat=cvfit$lambda.min
  fit=glmnet(x,y,standardize=standardize)
  yhat=predict(fit,x,s=lamhat)
  nz=sum(predict(fit,s=lamhat, type="coef")!=0)
  sigma=sqrt(sum((y-yhat)^2)/(length(y)-nz-1))
  return(list(sigmahat=sigma, df=nz))
}

# Update the QR factorization, after a column has been
# added. Here Q1 is m x n, Q2 is m x k, and R is n x n.

updateQR <- function(Q1,Q2,R,col) {
  m = nrow(Q1)
  n = ncol(Q1)
  k = ncol(Q2)
  
  a = .C("update1",
         Q2=as.double(Q2),
         w=as.double(t(Q2)%*%col),
         m=as.integer(m),
         k=as.integer(k),
         dup=FALSE,
         package="selectiveInference")
  
  Q2 = matrix(a$Q2,nrow=m)
  w = c(t(Q1)%*%col,a$w)
  
  # Re-structure: delete a column from Q2, add one to
  # Q1, and expand R
  Q1 = cbind(Q1,Q2[,1])
  Q2 = Q2[,-1,drop=FALSE]
  R = rbind(R,rep(0,n))
  R = cbind(R,w[Seq(1,n+1)])
  
  return(list(Q1=Q1,Q2=Q2,R=R))
}

# Moore-Penrose pseudo inverse for symmetric matrices

pinv <- function(A, tol=.Machine$double.eps) {
  e = eigen(A)
  v = Re(e$vec)
  d = Re(e$val)
  d[d > tol] = 1/d[d > tol]
  d[d < tol] = 0
  if (length(d)==1) return(v*d*v)
  else return(v %*% diag(d) %*% t(v))
}

##############################

# Assuming that grid is in sorted order from smallest to largest,
# and vals are monotonically increasing function values over the
# grid, returns the grid end points such that the corresponding
# vals are approximately equal to {val1, val2}

grid.search <- function(grid, fun, val1, val2, gridpts=100, griddepth=2) {
  n = length(grid)
  vals = fun(grid)
  
  ii = which(vals >= val1)
  jj = which(vals <= val2)
  if (length(ii)==0) return(c(grid[n],Inf))   # All vals < val1
  if (length(jj)==0) return(c(-Inf,grid[1]))  # All vals > val2
  # RJT: the above logic is correct ... but for simplicity, instead,
  # we could just return c(-Inf,Inf) 
  
  i1 = min(ii); i2 = max(jj)
  if (i1==1) lo = -Inf
  else lo = grid.bsearch(grid[i1-1],grid[i1],fun,val1,gridpts,
                         griddepth-1,below=TRUE)
  if (i2==n) hi = Inf
  else hi = grid.bsearch(grid[i2],grid[i2+1],fun,val2,gridpts,
                         griddepth-1,below=FALSE)
  return(c(lo,hi))
}

# Repeated bin search to find the point x in the interval [left, right]
# that satisfies f(x) approx equal to val. If below=TRUE, then we seek
# x such that the above holds and f(x) <= val; else we seek f(x) >= val.

grid.bsearch <- function(left, right, fun, val, gridpts=100, griddepth=1, below=TRUE) {
  n = gridpts
  depth = 1
  
  while (depth <= griddepth) {
    grid = seq(left,right,length=n)
    vals = fun(grid)
    
    if (below) {
      ii = which(vals >= val)
      if (length(ii)==0) return(grid[n])   # All vals < val (shouldn't happen)
      if ((i0=min(ii))==1) return(grid[1]) # All vals > val (shouldn't happen)
      left = grid[i0-1]
      right = grid[i0]
    }
    
    else {
      ii = which(vals <= val)
      if (length(ii)==0) return(grid[1])   # All vals > val (shouldn't happen)
      if ((i0=max(ii))==n) return(grid[n]) # All vals < val (shouldn't happen)
      left = grid[i0]
      right = grid[i0+1]
    }
    
    depth = depth+1
  }
  
  return(ifelse(below, left, right))
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector

tnorm.surv <- function(z, mean, sd, a, b, bits=NULL) {
  z = max(min(z,b),a)
  
  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1
  
  # Try the multi precision floating point calculation first
  o = is.finite(mean)
  mm = mean[o]
  pp = mpfr.tnorm.surv(z,mm,sd,a,b,bits) 
  
  # If there are any NAs, then settle for an approximation
  oo = is.na(pp)
  if (any(oo)) pp[oo] = bryc.tnorm.surv(z,mm[oo],sd,a,b)
  
  p[o] = pp
  return(p)
}

# Returns Prob(Z>z | Z in [a,b]), where mean cane be a vector, using
# multi precision floating point calculations thanks to the Rmpfr package

mpfr.tnorm.surv <- function(z, mean=0, sd=1, a, b, bits=NULL) {
  # If bits is not NULL, then we are supposed to be using Rmpf
  # (note that this was fail if Rmpfr is not installed; but
  # by the time this function is being executed, this should
  # have been properly checked at a higher level; and if Rmpfr
  # is not installed, bits would have been previously set to NULL)
  if (!is.null(bits)) {
    z = Rmpfr::mpfr((z-mean)/sd, precBits=bits)
    a = Rmpfr::mpfr((a-mean)/sd, precBits=bits)
    b = Rmpfr::mpfr((b-mean)/sd, precBits=bits)
    return(as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z))/
                        (Rmpfr::pnorm(b)-Rmpfr::pnorm(a))))
  }
  
  # Else, just use standard floating point calculations
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  return((pnorm(b)-pnorm(z))/(pnorm(b)-pnorm(a)))
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL, W Bryc
# Applied Mathematics and Computation
# Volume 127, Issues 23, 15 April 2002, Pages 365--374
# https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf

bryc.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  n = length(mean)
  
  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  p = (ff(z)-term2)/(term1-term2)
  
  # Sometimes the approximation can give wacky p-values,
  # outside of [0,1] ..
  #p[p<0 | p>1] = NA
  p = pmin(1,pmax(0,p))
  return(p)
}

ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
           (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
}

############## MODIFIED FUNCTIONS ###############

# Lasso inference function (for fixed lambda). Note: here we are providing inference
# for the solution of
# min 1/2 || y - \beta_0 - X \beta ||_2^2 + \lambda || \beta ||_1

oldFixedLassoInf <- function(x, y, beta, lambda, family=c("gaussian","binomial","cox"),intercept=TRUE, status=NULL,
                          sigma=NULL, alpha=0.1,
                          type=c("partial","full"), tol.beta=1e-5, tol.kkt=0.1,
                          gridrange=c(-100,100), bits=NULL, verbose=FALSE, offset_correction=TRUE) {
  
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
      hsigmaSinv <- pinv(hsigmaS) # solve(hsigmaS)
      
      # Approximate inverse covariance matrix for when (n < p) from lasso_Inference.R
      htheta <- InverseLinfty(hsigma, n, verbose=FALSE)
      
      # 0-padding matrix
      FS = rbind(diag(length(S)),matrix(0,pp-length(S),length(S)))
      ithetasigma = (diag(pp)-(htheta%*%hsigma))

      M <- (((htheta%*%t(Xordered))+ithetasigma%*%FS%*%hsigmaSinv%*%t(XS))/n)
      # vector which is offset for testing debiased beta's
      meanoffset <- -(((ithetasigma%*%FS%*%hsigmaSinv)%*%sign(hbetaS))*lambda/n)
      if (intercept == T) {
        M = M[-1,] # remove intercept row
        meanoffset = meanoffset[-1] # remove intercept element
      }
      if (offset_correction == FALSE) {
        meanoffset = 0 * meanoffset
      }
    } else if (type=="partial" || p > n) {
      xa = x[,vars,drop=F]
      M = pinv(crossprod(xa)) %*% t(xa)
      meanoffset = rep(0,k)
    } else {
      M = pinv(crossprod(x)) %*% t(x)
      M = M[vars,,drop=F]
      meanoffset = rep(0,k)
    }
    
    for (j in 1:k) {
      if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))
      
      vj = M[j,]
      mj = sqrt(sum(vj^2))
      vj = vj / mj        # Standardize (divide by norm of vj)
      sign[j] = sign(sum(vj*y))
      vj = sign[j] * vj

      a = poly.pval(y,G,u,vj,offset=meanoffset[j],sigma,bits)
      pv[j] = a$pv
      vlo[j] = a$vlo * mj # Unstandardize (mult by norm of vj)
      vup[j] = a$vup * mj # Unstandardize (mult by norm of vj)
      vmat[j,] = vj * mj * sign[j]  # Unstandardize (mult by norm of vj)
      
      a = poly.int(y,G,u,vj,offset=meanoffset[j],sigma,alpha,gridrange=gridrange,
                   flip=(sign[j]==-1),bits=bits)
      ci[j,] = a$int * mj # Unstandardize (mult by norm of vj)
      tailarea[j,] = a$tailarea
    }
    
    out = list(type=type,lambda=lambda,pv=pv,ci=ci,
               tailarea=tailarea,vlo=vlo,vup=vup,vmat=vmat,y=y,
               vars=vars,sign=sign,sigma=sigma,alpha=alpha,
               sd=sigma*sqrt(rowSums(vmat^2)),
               coef0=vmat%*%y,
               call=this.call,M=M)
    class(out) = "fixedLassoInf"
    return(out)
  }
}


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


# Main p-value function

poly.pval <- function(y, G, u, v, sigma, offset=0, bits=NULL) {
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))
  
  pv = tnorm.surv(z,0-offset,sd,vlo,vup,bits)
  return(list(pv=pv,vlo=vlo,vup=vup))
}

# Main confidence interval function

poly.int <- function(y, G, u, v, sigma, alpha, offset=0, gridrange=c(-100,100),
                     gridpts=100, griddepth=2, flip=FALSE, bits=NULL) {
  
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))
  
  xg = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
  fun = function(x) { tnorm.surv(z,x-offset,sd,vlo,vup,bits) }
  
  int = grid.search(xg,fun,alpha/2,1-alpha/2,gridpts,griddepth)
  tailarea = c(fun(int[1]),1-fun(int[2]))
  
  if (flip) {
    int = -int[2:1]
    tailarea = tailarea[2:1]
  }
  
  return(list(int=int,tailarea=tailarea))
}

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

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

checkargs.xy <- function(x, y) {
  if (missing(x)) stop("x is missing.")
  if (is.null(x) || !is.matrix(x)) stop("x must be a matrix.")
  if (missing(y)) stop("y is missing.")
  if (is.null(y) || !is.numeric(y)) stop("y must be numeric.")
  if (ncol(x) == 0) stop("There must be at least one predictor [must have ncol(x) > 0].")
  if (checkcols(x)) stop("x cannot have duplicate columns.")
  if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 0].")
  if (length(y)!=nrow(x)) stop("Dimensions don't match [length(y) != nrow(x)].")
}

checkargs.misc <- function(sigma=NULL, alpha=NULL, lambda=NULL, k=NULL) {
  if (!is.null(sigma) && sigma<=0) stop("sigma must be > 0")
  if (!is.null(lambda) && lambda<0) stop("lambda must be >= 0")
  if (!is.null(alpha) && (alpha <=0 || alpha >= 1)) stop("alpha must be between 0 and 1")
  if (!is.null(k) && (k < 1 || k!=floor(k))) stop("k must be an integer >= 1")
}

# Make sure that no two columms of A are the same
# (this works with probability one).

checkcols <- function(A) {
  b = rnorm(nrow(A))
  a = sort(t(A)%*%b)
  if (any(diff(a)==0)) return(TRUE)
  return(FALSE)
}

checkargs2= function(y,x=NULL,bhat=NULL,sigma=NULL,lambda=0,alpha=.1,tol.beta=1e-5,tol.kkt=0.01,nsteps=NULL
,larfit=NULL, fsfit=NULL, k=NULL,bh.q=NULL){
  
    # checks arguments for all user-callable functions.
  
 if(!is.null(x)){
    if(!is.matrix(x)) stop("x must be matrix")
    if(ncol(x)<2) stop("x must have at least 2 columns")
    if(length(y)!=nrow(x)) stop("Number of rows of x must equal length of y")
}
 
    if(!is.null(bhat)){
         if(length(bhat)!=ncol(x))  stop("bhat must have length equal to the number of cols of x")
     }
    if(!is.null(sigma)){
         if(sigma<=0) stop("sigma must be gt 0")
     }
    if(lambda<0) stop("lambda must be non-negative")
    if(alpha<=0 | alpha >=1) stop("alpha must be between 0 and 1")
      if(tol.beta<=0) stop("tol.beta must be gt 0")
    if(tol.kkt<=0) stop("tol.kkt must be gt 0")
    if(!is.null(nsteps)){
        if(nsteps<1 | nsteps> min(nrow(x),ncol(x))) stop("nsteps must be between 1 and min(nrow(x),ncol(x))")
        if(!is.wholenumber(nsteps)) stop("nsteps must be an integer")
    }
  
       if(!is.null(larfit)){
             if(class(larfit)!="lar") stop("larfit must be a lar object returned by lar (not lars)")
         }
    if(!is.null(fsfit)){
             if(class(fsfit)!="forwardStep") stop("fsfit must be an object returned by forwardStep")
         }
if(!is.null(bh.q)){
  if(bh.q<=0 | bh.q >=1) stop("bh.q must be between 0 and 1")
   }
   if(!is.null(k)){
        if(k<1 | k >length(y))  stop("k must be between 1 and length(y)")
    } 
}

estimateSigma=function(x,y){
    if(nrow(x)<10) stop("Number of observations must be at least 10 to run estimateSigma")
    cvfit=cv.glmnet(x,y,standardize=F)
    lamhat=cvfit$lambda.min
    fit=glmnet(x,y,standardize=F)
    yhat=predict(fit,x,s=lamhat)
    nz=sum(predict(fit,s=lamhat, type="coef")!=0)
    sigma=sqrt(sum((y-yhat)^2)/(length(y)-nz-1))
    return(list(sigmahat=sigma, df=nz))
}


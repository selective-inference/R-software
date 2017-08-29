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

coef.interpolate <- function(beta, s, knots, decreasing=TRUE) {
  # Sort the s values
  o = order(s,decreasing=decreasing)
  s = s[o]

  k = length(s)
  mat = matrix(rep(knots,each=k),nrow=k)
  if (decreasing) b = s >= mat
  else b = s <= mat
  blo = max.col(b,ties.method="first")
  bhi = pmax(blo-1,1)

  i = bhi==blo
  p = numeric(k)
  p[i] = 0
  p[!i] = ((s-knots[blo])/(knots[bhi]-knots[blo]))[!i]

  beta = t((1-p)*t(beta[,blo,drop=FALSE]) + p*t(beta[,bhi,drop=FALSE]))
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
  if (!is.null(lambda) && length(lambda) != 1) stop("lambda must be a single number")
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

  a = update1_(as.matrix(Q2), t(Q2)%*%col, m, k) # Rcpp call

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

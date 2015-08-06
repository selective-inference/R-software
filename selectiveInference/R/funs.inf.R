# Main p value function

Gv.pval <- function(y, G, v, sigma) {
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  r = G %*% v
  vec = -G %*% (y - v*sum(v*y)/vv) / r * vv
  vlo = suppressWarnings(max(vec[r>0]))
  vup = suppressWarnings(min(vec[r<0]))
  
  pv = tnorm.surv(z,0,sd,vlo,vup)
  return(list(pv=pv,vlo=vlo,vup=vup))
}

# Main confidence interval function

Gv.int <- function(y, G, v, sigma, alpha, gridfac=25, gridpts=1000,
                   flip=FALSE) {
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  r = G %*% v
  vec = -G %*% (y - v*sum(v*y)/vv) / r * vv
  vlo = suppressWarnings(max(vec[r>0]))
  vup = suppressWarnings(min(vec[r<0]))

  fun = function(mu) return(tnorm.surv(z,mu,sd,vlo,vup))
  xl = z-gridfac*sd
  xr = z+gridfac*sd
  lo = grid.search(fun,alpha/2,xl,xr,gridpts,left.side=TRUE)
  hi = grid.search(fun,1-alpha/2,xl,xr,gridpts,left.side=FALSE)
  int = c(lo,hi)
  tailarea = c(fun(lo),1-fun(hi)) 

  if (flip) {
    int = -int[2:1]
    tailarea = tailarea[2:1]
  }
  
  return(list(int=int,tailarea=tailarea))
}

##############################

grid.search <- function(fun, val, xl, xr, gridpts=1000, left.side=TRUE) {
  xg = seq(xl,xr,length=gridpts)
  vals = fun(xg)
  # return(xg[which.min(abs(vals-val))])
  
  if (left.side) {
    ii = which(vals<=val)
    if (length(ii)==0) return(-Inf)
    else return(xg[max(ii)])
  }
  else {
    ii = which(vals>=val)
    if (length(ii)==0) return(Inf)
    else return(xg[min(ii)])
  }
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector

tnorm.surv <- function(z,mean,sd,a,b) {
  z = max(min(z,b),a)

  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1
  
  o = is.finite(mean)
  p[o] = rob.tnorm.surv(z,mean[o],sd,a,b)  
  return(p)
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector,
# based on
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL, W Bryc
# Applied Mathematics and Computation
# Volume 127, Issues 2–3, 15 April 2002, Pages 365–374
# https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf

rob.tnorm.surv <- function(z, mean=0, sd=1, a, b){
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd

  term1 = exp(z*z)
  o = a >- Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = 0
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  p = (ff(z)-term2)/(term1-term2)
  return(pmax(pmin(p,1),0))
}

rob.tnorm.surv.old <- function(z ,mean=0, sd=1, a, b){
  # First try standard formula
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  p = (pnorm(b)-pnorm(z))/(pnorm(b)-pnorm(a))
  o = is.na(p)
  
  # If there are NAs, apply modified method
  if (any(o)) {
    zz = z[o]
    aa = a[o]
    bb = b[o]
    term1 = exp(zz*zz)
    oo = aa >- Inf
    term1[oo] = ff(aa[oo])*exp(-(aa[oo]^2-zz[oo]^2)/2)
    term2 = 0
    ooo = bb < Inf
    term2[ooo] = ff(bb[ooo])*exp(-(bb[ooo]^2-zz[ooo]^2)/2)
    p[o] = (ff(zz)-term2)/(term1-term2)
  }
  
  return(p)
}

ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
         (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
}

##############################

forwardStop <- function(pv, alpha=.10){
  if (alpha<0 || alpha>1) stop("alpha must be in [0,1]")
  if (min(pv,na.rm=T)<0 || max(pv,na.rm=T)>1) stop("pvalues must be in [0,1]")
  val=-(1/(1:length(pv)))*cumsum(log(1-pv))
  oo = which(val <= alpha)
  if (length(oo)==0) out=0
  else out = oo[length(oo)]
  return(out)
}


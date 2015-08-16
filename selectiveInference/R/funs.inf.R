# Main p-value function

poly.pval <- function(y, G, u, v, sigma) {
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))

  pv = tnorm.surv(z,0,sd,vlo,vup)
  return(list(pv=pv,vlo=vlo,vup=vup))
}

# Main confidence interval function

poly.int <- function(y, G, u, v, sigma, alpha, gridrange=c(-100,100),
                     gridpts=1000, flip=FALSE) {

  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))

  xg = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
  vals = tnorm.surv(z,xg,sd,vlo,vup)
  int = grid.search(xg,vals,alpha/2,1-alpha/2)
  tailarea = c(tnorm.surv(z,int[1],sd,vlo,vup),
    1-tnorm.surv(z,int[2],sd,vlo,vup))

  if (flip) {
    int = -int[2:1]
    tailarea = tailarea[2:1]
  }
  
  return(list(int=int,tailarea=tailarea))
}

##############################

# Assuming that grid is in sorted order from smallest to largest,
# and vals are monotonically increasing function values over the
# grid, returns the grid end points such that the corresponding
# vals are approximately equal to {left, right}

grid.search <- function(grid, vals, left, right) {
  n = length(grid)
  il = which(vals >= left)
  ir = which(vals <= right)
  if (length(il)==0) return(c(grid[n],Inf))   # All vals < left
  if (length(ir)==0) return(c(-Inf,grid[1]))  # All vals > right
  ## RJT QUESTION: should we just return c(-Inf,Inf) for the
  ## above two cases?? The above logic seems correct ...

  i1 = min(il); i2 = max(ir)
  if (i1==1) lo = -Inf
  else lo = grid[i1-1]
  if (i2==n) hi = Inf
  else hi = grid[i2+1]
  return(c(lo,hi))
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector

tnorm.surv <- function(z, mean, sd, a, b) {
  z = max(min(z,b),a)
  
  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1
  
  o = is.finite(mean)
  p[o] = bryc.tnorm.surv(z,mean[o],sd,a,b)
  #p[o] = gsell.tnorm.surv(z,mean[o],sd,a,b)
  return(p)
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

# Return Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# Riemann approximation tricks, by Max G'Sell

gsell.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  return(max.approx.frac(a/sd,b/sd,z/sd,mean/sd))
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

##############################

aicStop <- function(x, y, action, df, sigma, mult=2, ntimes=2) {
  n = length(y)
  k = length(action)
  aic = numeric(k)
  G = matrix(0,nrow=0,ncol=n)
  u = numeric(0)
  count = 0
  
  for (i in 1:k) {
    A = action[1:i]
    aic[i] = sum(lsfit(x[,A],y,intercept=F)$res^2) + mult*sigma^2*df[i]

    j = action[i]
    if (i==1) xtil = x[,j]
    else xtil = lsfit(x[,action[1:(i-1)]],x[,j],intercept=F)$res
    s = sign(sum(xtil*y))
    
    if (i==1 || aic[i] <= aic[i-1]) {
      G = rbind(G,s*xtil/sqrt(sum(xtil^2)))
      u = c(u,sqrt(mult)*sigma)
      count = 0
    }

    else {
      G = rbind(G,-s*xtil/sqrt(sum(xtil^2)))
      u = c(u,-sqrt(mult)*sigma)
      count = count+1
      if (count == ntimes) break
    }
  }

  if (i < k) {
    khat = i - ntimes
    aic = aic[1:i]
  }
  else khat = k
  return(list(khat=khat,G=G,u=u,aic=aic))
}

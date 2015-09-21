# We compute the forward stepwise regression (FS) path given 
# a response vector y and predictor matrix x.  We assume
# that x has columns in general position.

fs <- function(x, y, maxsteps=2000, intercept=TRUE, normalize=TRUE,
               verbose=FALSE) {

  this.call = match.call()
  checkargs.xy(x=x,y=y)
  
  # Center and scale, etc.
  obj = standardize(x,y,intercept,normalize)
  x = obj$x
  y = obj$y
  bx = obj$bx
  by = obj$by
  sx = obj$sx
  n = nrow(x)
  p = ncol(x)

  #####
  # To keep consistent with the lar function, we parametrize
  # so that the first step has all zero coefficients,
  # Also, an interesting note: the effective "lambda" (maximal
  # correlation with the residual) may increase with stepwise!
  # So we don't keep track of it

  #####
  # Find the first variable to enter and its sign
  xx = scale(x,center=F,scale=sqrt(colSums(x^2)))
  uhat = t(xx)%*%y
  ihit = which.max(abs(uhat))   # Hitting coordinate
  s = Sign(uhat[ihit])          # Sign

  if (verbose) {
    cat(sprintf("1. Adding variable %i, |A|=%i...",ihit,1))
  }
  
  # Now iteratively find the new FS estimates

  # Things to keep track of, and return at the end
  buf = min(maxsteps,500)
  action = numeric(buf)      # Actions taken
  df = numeric(buf)          # Degrees of freedom
  beta = matrix(0,p,buf)     # FS estimates
  
  action[1] = ihit
  df[1] = 0
  beta[,1] = 0

  # Gamma matrix!
  gbuf = max(2*p*3,2000)     # Space for 3 steps, at least
  gi = 0
  Gamma = matrix(0,gbuf,n)
  Gamma[gi+Seq(1,p-1),] = t(s*xx[,ihit]+xx[,-ihit]); gi = gi+p-1
  Gamma[gi+Seq(1,p-1),] = t(s*xx[,ihit]-xx[,-ihit]); gi = gi+p-1
  Gamma[gi+1,] = t(s*xx[,ihit]); gi = gi+1

  # nk
  nk = numeric(buf)
  vreg = matrix(0,buf,n)
  nk[1] = gi
  vreg[1,] = s*x[,ihit] / sum(x[,ihit]^2)

  # Other things to keep track of, but not return
  r = 1                      # Size of active set
  A = ihit                   # Active set
  I = Seq(1,p)[-ihit]        # Inactive set
  X1 = x[,ihit,drop=FALSE]   # Matrix X[,A]
  X2 = x[,-ihit,drop=FALSE]  # Matrix X[,I]
  k = 2                      # What step are we at?

  # Compute a skinny QR decomposition of X1
  obj = qr(X1)
  Q = qr.Q(obj,complete=TRUE)
  Q1 = Q[,1,drop=FALSE];
  Q2 = Q[,-1,drop=FALSE]
  R = qr.R(obj)
  
  # Throughout the algorithm, we will maintain
  # the decomposition X1 = Q1*R. Dimenisons:
  # X1: n x r
  # Q1: n x r
  # Q2: n x (n-r)
  # R:  r x r
    
  while (k<=maxsteps) {
    ##########
    # Check if we've reached the end of the buffer
    if (k > length(action)) {
      buf = length(action)
      action = c(action,numeric(buf))
      df = c(df,numeric(buf))
      beta = cbind(beta,matrix(0,p,buf))
      nk = c(nk,numeric(buf))
      vreg = rbind(vreg,matrix(0,buf,n))
    }

    # Key quantities for the next entry
    a = backsolve(R,t(Q1)%*%y)
    X2perp = X2 - X1 %*% backsolve(R,t(Q1)%*%X2)
    xx = scale(X2perp,center=F,scale=sqrt(colSums(X2perp^2)))
    aa = as.numeric(t(xx)%*%y)
    
    # If the inactive set is empty, nothing will hit
    if (r==min(n-intercept,p)) break

    # Otherwise find the next hitting time
    else {
      shits = Sign(aa)
      hits = shits * aa
      ihit = which.max(hits)
      shit = shits[ihit]
    }
    
    # Record the solution
    action[k] = I[ihit]
    df[k] = r
    beta[A,k] = a
        
    # Gamma matrix!
    if (gi + 2*p > nrow(Gamma)) Gamma = rbind(Gamma,matrix(0,2*p+gbuf,n))
    xx = t(shits*t(xx))
    Gamma[gi+Seq(1,p-r),] = t(xx); gi = gi+p-r
    Gamma[gi+Seq(1,p-r-1),] = t(xx[,ihit]-xx[,-ihit]); gi = gi+p-r-1
    Gamma[gi+1,] = t(xx[,ihit]); gi = gi+1

    # nk, regression contrast
    nk[k] = gi
    vreg[k,] = shit*X2perp[,ihit] / sum(X2perp[,ihit]^2)

    # Update all of the variables
    r = r+1
    A = c(A,I[ihit])
    I = I[-ihit]
    s = c(s,shit)
    X1 = cbind(X1,X2[,ihit])
    X2 = X2[,-ihit,drop=FALSE]

    # Update the QR decomposition
    obj = updateQR(Q1,Q2,R,X1[,r])
    Q1 = obj$Q1
    Q2 = obj$Q2
    R = obj$R
     
    if (verbose) {
      cat(sprintf("\n%i. Adding variable %i, |A|=%i...",k,A[r],r))
    }
            
    # Step counter
    k = k+1
  }

  # Trim
  action = action[Seq(1,k-1)]
  df = df[Seq(1,k-1),drop=FALSE]
  beta = beta[,Seq(1,k-1),drop=FALSE]
  Gamma = Gamma[Seq(1,gi),,drop=FALSE]
  nk = nk[Seq(1,k-1)]
  vreg = vreg[Seq(1,k-1),,drop=FALSE]
  
  # If we reached the maximum number of steps
  if (k>maxsteps) {
    if (verbose) {
      cat(sprintf("\nReached the maximum number of steps (%i),",maxsteps))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
    bls = NULL
  }

  # Otherwise, note that we completed the path
  else {
    completepath = TRUE
    
    # Record the least squares solution. Note that
    # we have already computed this
    bls = rep(0,p)
    bls[A] = a
  }

  if (verbose) cat("\n")
  
  # Adjust for the effect of centering and scaling
  if (intercept) df = df+1
  if (normalize) beta = beta/sx
  if (normalize && completepath) bls = bls/sx
  
  # Assign column names
  colnames(beta) = as.character(Seq(1,k-1))

  out = list(action=action,sign=s,df=df,beta=beta,
    completepath=completepath,bls=bls,
    Gamma=Gamma,nk=nk,vreg=vreg,x=x,y=y,bx=bx,by=by,sx=sx,
    intercept=intercept,normalize=normalize,call=this.call) 
  class(out) = "fs"
  return(out)
}

##############################

# Coefficient function for fs

coef.fs <- function(object, s, ...) {
  if (object$completepath) {
    k = length(object$action)+1
    beta = cbind(object$beta,object$bls)
  } else {
    k = length(object$action)
    beta = object$beta
  }
  
  if (min(s)<0 || max(s)>k) stop(sprintf("s must be between 0 and %i",k))
  knots = 1:k
  dec = FALSE
  return(coef.interpolate(beta,s,knots,dec))
}

# Prediction function for fs

predict.fs <- function(object, newx, s, ...) {
  beta = coef.fs(object,s)
  if (missing(newx)) newx = scale(object$x,FALSE,1/object$sx)
  else newx = scale(newx,object$bx,FALSE)
  return(newx %*% beta + object$by)
}

##############################

# FS inference function

fsInf <- function(obj, sigma=NULL, alpha=0.1, k=NULL, type=c("active","all","aic"), 
                  gridrange=c(-100,100), bits=NULL, mult=2, ntimes=2, verbose=FALSE) {
  
  this.call = match.call()
  type = match.arg(type)
  checkargs.misc(sigma=sigma,alpha=alpha,k=k,
                 gridrange=gridrange,mult=mult,ntimes=ntimes)
  if (class(obj) != "fs") stop("obj must be an object of class fs")
  if (is.null(k) && type=="active") k = length(obj$action)
  if (is.null(k) && type=="all") stop("k must be specified when type = all")
  if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
    warning("Package Rmpfr is not installed, reverting to standard precision")
    bits = NULL
  }
  
  k = min(k,length(obj$action)) # Round to last step
  x = obj$x
  y = obj$y
  p = ncol(x)
  n = nrow(x)
  G = obj$Gamma
  nk = obj$nk
  sx = obj$sx

  if (is.null(sigma)) {
    if (n >= 2*p) {
      oo = obj$intercept
      sigma = sqrt(sum(lsfit(x,y,intercept=oo)$res^2)/(n-p-oo))
    }
    else {
      sigma = sd(y)
      warning(paste(sprintf("p > n/2, and sd(y) = %0.3f used as an estimate of sigma;",sigma),
                    "you may want to use the estimateSigma function"))
    }
  }

  khat = NULL
  
  if (type == "active") {
    pv = vlo = vup = numeric(k) 
    vmat = matrix(0,k,n)
    ci = tailarea = matrix(0,k,2)
    vreg = obj$vreg[1:k,,drop=FALSE]
    sign = obj$sign[1:k]
    vars = obj$action[1:k]

    for (j in 1:k) {
      if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))

      Gj = G[1:nk[j],]
      uj = rep(0,nk[j])
      vj = vreg[j,]
      mj = sqrt(sum(vj^2)) 
      vj = vj / mj              # Standardize (divide by norm of vj)
      a = poly.pval(y,Gj,uj,vj,sigma,bits)
      pv[j] = a$pv
      sxj = sx[vars[j]]
      vlo[j] = a$vlo * mj / sxj # Unstandardize (mult by norm of vj / sxj)
      vup[j] = a$vup * mj / sxj # Unstandardize (mult by norm of vj / sxj)
      vmat[j,] = vj * mj / sxj  # Unstandardize (mult by norm of vj / sxj)
  
      a = poly.int(y,Gj,uj,vj,sigma,alpha,gridrange=gridrange,
        flip=(sign[j]==-1),bits=bits)
      ci[j,] = a$int * mj / sxj # Unstandardize (mult by norm of vj / sxj)
      tailarea[j,] = a$tailarea
    }

    khat = forwardStop(pv,alpha)
  }
  
  else {
    if (type == "aic") {
      out = aicStop(x,y,obj$action[1:k],obj$df[1:k],sigma,mult,ntimes)
      khat = out$khat
      m = out$stopped * ntimes
      G = rbind(out$G,G[1:nk[khat+m],])  # Take ntimes more steps past khat
      u = c(out$u,rep(0,nk[khat+m]))     # (if we need to)
      kk = khat
    }
    else {
      G = G[1:nk[k],]
      u = rep(0,nk[k])
      kk = k
    }
    
    pv = vlo = vup = numeric(kk) 
    vmat = matrix(0,kk,n)
    ci = tailarea = matrix(0,kk,2)
    sign = numeric(kk)
    vars = obj$action[1:kk]
    xa = x[,vars]
    M = pinv(crossprod(xa)) %*% t(xa)
    
    for (j in 1:kk) {
      if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))

      vj = M[j,]
      mj = sqrt(sum(vj^2))
      vj = vj / mj              # Standardize (divide by norm of vj)
      sign[j] = sign(sum(vj*y))
      vj = sign[j] * vj
      Gj = rbind(G,vj)
      uj = c(u,0)

      a = poly.pval(y,Gj,uj,vj,sigma,bits)
      pv[j] = a$pv
      sxj = sx[vars[j]]
      vlo[j] = a$vlo * mj / sxj # Unstandardize (mult by norm of vj / sxj)
      vup[j] = a$vup * mj / sxj # Unstandardize (mult by norm of vj / sxj)
      vmat[j,] = vj * mj / sxj  # Unstandardize (mult by norm of vj / sxj)

      a = poly.int(y,Gj,uj,vj,sigma,alpha,gridrange=gridrange,
        flip=(sign[j]==-1),bits=bits)
      ci[j,] = a$int * mj / sxj # Unstandardize (mult by norm of vj / sxj)
      tailarea[j,] = a$tailarea
    }
  }
  
  out = list(type=type,k=k,khat=khat,pv=pv,ci=ci,
    tailarea=tailarea,vlo=vlo,vup=vup,vmat=vmat,y=y,
    vars=vars,sign=sign,sigma=sigma,alpha=alpha,
    call=this.call)
  class(out) = "fsInf"
  return(out)
}

##############################


##############################

print.fs <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call)
  
  cat("\nSequence of FS moves:\n")
  nsteps = length(x$action)
  tab = cbind(1:nsteps,x$action,x$sign)
  colnames(tab) = c("Step","Var","Sign")
  rownames(tab) = rep("",nrow(tab))
  print(tab)
  invisible()
}

print.fsInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))

  if (x$type == "active") {
    cat(sprintf("\nSequential testing results with alpha = %0.3f\n",x$alpha))
    tab = cbind(1:length(x$pv),x$vars,
      round(x$sign*x$vmat%*%x$y,3),
      round(x$sign*x$vmat%*%x$y/(x$sigma*sqrt(rowSums(x$vmat^2))),3),
      round(x$pv,3),round(x$ci,3))
    colnames(tab) = c("Step", "Var", "Coef", "Z-score", "P-value",
              "LowConfPt", "UpConfPt")
    if (tailarea) {
      tab = cbind(tab,round(x$tailarea,3))
      colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
    }
    rownames(tab) = rep("",nrow(tab))
    print(tab)

    cat(sprintf("\nEstimated stopping point from ForwardStop rule = %i\n",x$khat))
  }

  else if (x$type == "all") {
    cat(sprintf("\nTesting results at step = %i, with alpha = %0.3f\n",x$k,x$alpha))
    tab = cbind(x$vars,
      round(x$sign*x$vmat%*%x$y,3),
      round(x$sign*x$vmat%*%x$y/(x$sigma*sqrt(rowSums(x$vmat^2))),3),
      round(x$pv,3),round(x$ci,3))
    colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
    if (tailarea) {
      tab = cbind(tab,round(x$tailarea,3))
      colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
    }
    rownames(tab) = rep("",nrow(tab))
    print(tab)
  }

  else if (x$type == "aic") {
    cat(sprintf("\nTesting results at step = %i, with alpha = %0.3f\n",x$khat,x$alpha))
    tab = cbind(x$vars,
      round(x$sign*x$vmat%*%x$y,3),
      round(x$sign*x$vmat%*%x$y/(x$sigma*sqrt(rowSums(x$vmat^2))),3),
      round(x$pv,3),round(x$ci,3))
    colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
    if (tailarea) {
      tab = cbind(tab,round(x$tailarea,3))
      colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
    }
    rownames(tab) = rep("",nrow(tab))
    print(tab)
    
    cat(sprintf("\nEstimated stopping point from AIC rule = %i\n",x$khat))
  }

  invisible()
}

plot.fs <- function(x, breaks=TRUE, omit.zeros=TRUE, var.labels=TRUE, ...) {
  if (x$completepath) {
    k = length(x$action)+1
    beta = cbind(x$beta,x$bls)
  } else {
    k = length(x$action)
    beta = x$beta
  }
  p = nrow(beta)
  
  xx = 1:k
  xlab = "Step"
  
  if (omit.zeros) {
    inds = matrix(FALSE,p,k)
    for (i in 1:k) {
      inds[i,] = beta[i,]!=0 | c(diff(beta[i,]!=0),F) | c(F,diff(beta[i,]!=0))
    }
    beta[!inds] = NA
  }

  plot(c(),c(),xlim=range(xx,na.rm=T),ylim=range(beta,na.rm=T),
       xlab=xlab,ylab="Coefficients",main="Forward stepwise path",...)
  abline(h=0,lwd=2)
  matplot(xx,t(beta),type="l",lty=1,add=TRUE)
  if (breaks) abline(v=x,lty=2)
  if (var.labels) axis(4,at=beta[,k],labels=1:p,cex=0.8,adj=0) 
  invisible()
}


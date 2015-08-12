# We compute the forward stepwise regression (FS) path given 
# a response vector y and predictor matrix x.  We assume
# that x has columns in general position.

fs <- function(x, y, maxsteps=2000, verbose=FALSE, intercept=TRUE,
               normalize=TRUE) {

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
  Gamma = matrix(0,0,n)
  if (p>1) Gamma = rbind(Gamma,t(s*xx[,ihit]+xx[,-ihit]),t(s*xx[,ihit]-xx[,-ihit]))
  Gamma = rbind(Gamma,t(s*xx[,ihit]))
  nk = nrow(Gamma)

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
    }

    # Key quantities for the next entry
    a = backsolve(R,t(Q1)%*%y)
    mat = X2 - X1 %*% backsolve(R,t(Q1)%*%X2)
    xx = scale(mat,center=F,scale=sqrt(colSums(mat^2)))
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
    xx = t(shits*t(xx))
    Gamma = rbind(Gamma,t(xx))
    if (ncol(xx)>1) Gamma = rbind(Gamma,t(xx[,ihit]-xx[,-ihit]))
    Gamma = rbind(Gamma,t(xx[,ihit]))
    nk = c(nk,nrow(Gamma))

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
    Gamma=Gamma,nk=nk,x=x,y=y,bx=bx,by=by,sx=sx,
    intercept=intercept,normalize=normalize,call=this.call) 
  class(out) = "fs"
  return(out)
}

##############################

# Coefficient function for fs

coef.fs <- function(obj, s) {
  if (obj$completepath) {
    k = length(obj$action)+1
    beta = cbind(obj$beta,obj$bls)
  } else {
    k = length(obj$action)
    beta = obj$beta
  }
  
  if (min(s)<0 || max(s)>k) stop(sprintf("s must be between 0 and %i",k))
  knots = 1:k
  dec = FALSE
  return(coef.interpolate(beta,s,knots,dec))
}

# Prediction function for fs

predict.fs <- function(obj, newx, s) {
  beta = coef.fs(obj,s)
  if (missing(newx)) newx = scale(obj$x,FALSE,1/obj$sx)
  else newx = scale(newx,obj$bx,FALSE)
  return(newx %*% beta + obj$by)
}

##############################

# FS inference function

fsInf <- function(obj, sigma=NULL, alpha=0.1, k=NULL, type=c("active","all","aic"), 
                  gridfac=25, gridpts=1000, mult=2, ntimes=2) {
  
  this.call = match.call()
  type = match.arg(type)
  checkargs.misc(sigma=sigma,alpha=alpha,k=k)
  if (class(obj) != "fs") stop("obj must be an object of class fs")
  if (is.null(k) && type=="active") k = length(obj$action)
  if (is.null(k) && type=="all") stop("k must be specified when type = all")
  
  k = min(k,length(obj$action)) # Round to last step
  x = obj$x
  y = obj$y
  p = ncol(x)
  n = nrow(x)
  G = obj$Gamma
  nk = obj$nk

  if (is.null(sigma)) {
    if (n < 2*p) sigma = sd(y)
    else sigma = sqrt(sum(lsfit(x,y,intercept=F)$res^2)/(n-p))
  }

  khat = NULL
  
  if (type == "active") {
    pv = vlo = vup = numeric(k) 
    vmat = matrix(0,k,n)
    ci = tailarea = matrix(0,k,2)
    sign = obj$sign[1:k]
    vars = obj$action[1:k]

    for (j in 1:k) {
      Gj = G[1:nk[j],]
      uj = rep(0,nk[j])
      vj = G[nk[j],]
      vj = vj / sqrt(sum(vj^2))
      a = poly.pval(y,Gj,uj,vj,sigma)
      pv[j] = a$pv
      vlo[j] = a$vlo
      vup[j] = a$vup
      vmat[j,] = vj
    
      a = poly.int(y,Gj,uj,vj,sigma,alpha,gridfac=gridfac,gridpts=gridpts,
        flip=(sign[j]==-1))
      ci[j,] = a$int
      tailarea[j,] = a$tailarea
    }

    khat = forwardStop(pv,alpha)
  }
  
  else {
    if (type == "aic") {
      out = aicStop(x,y,obj$action[1:k],obj$df[1:k],sigma,mult,ntimes)
      khat = out$khat
      GG = out$G
      uu = out$u
      kk = khat
    }
    else {
      GG = matrix(0,0,n)
      uu = c()
      kk = k
    }
    
    pv = vlo = vup = numeric(kk) 
    vmat = matrix(0,kk,n)
    ci = tailarea = matrix(0,kk,2)
    sign = numeric(kk)
    vars = obj$action[1:kk]

    G = rbind(GG,G[1:nk[kk],])
    u = c(uu,rep(0,nk[kk]))
    xa = x[,vars]
    M = solve(crossprod(xa),t(xa))
    
    for (j in 1:kk) {
      vj = M[j,]
      sign[j] = sign(sum(vj*y))
      
      vj = vj / sqrt(sum(vj^2))
      vj = sign[j] * vj
      Gj = rbind(G,vj)
      uj = c(u,0)

      a = poly.pval(y,Gj,uj,vj,sigma)
      pv[j] = a$pv
      vlo[j] = a$vlo
      vup[j] = a$vup
      vmat[j,] = vj

      a = poly.int(y,Gj,uj,vj,sigma,alpha,gridfac=gridfac,gridpts=gridpts,
        flip=(sign[j]==-1))
      ci[j,] = a$int
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

print.fsInf <- function(x) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))

  if (x$type == "active") {
    cat(sprintf("\nSequential testing results with alpha = %0.3f\n",x$alpha))
    tab = cbind(1:length(x$pv),x$vars,
      round(x$sign*x$vmat%*%x$y,3),round(x$pv,3),round(x$ci,3),
      round(x$tailarea,3))
    colnames(tab) = c("Step", "Var", "Stdz Coef", "P-value", "LowConfPt",
              "UpConfPt", "LowArea", "UpArea")
    rownames(tab) = rep("",nrow(tab))
    print(tab)

    cat(sprintf("\nEstimated stopping point from ForwardStop rule = %i\n",x$khat))
  }

  else if (x$type == "all") {
    cat(sprintf("\nTesting results at step = %i, with alpha = %0.3f\n",x$k,x$alpha))
    tab = cbind(x$vars,round(x$sign*x$vmat%*%x$y,3),
      round(x$pv,3),round(x$ci,3),round(x$tailarea,3))
    colnames(tab) = c("Var", "Stdz Coef", "P-value", "LowConfPt", "UpConfPt",
              "LowArea", "UpArea")
    rownames(tab) = rep("",nrow(tab))
    print(tab)
  }

  else if (x$type == "aic") {
    cat(sprintf("\nTesting results at step = %i, with alpha = %0.3f\n",x$k,x$alpha))
    tab = cbind(x$vars,round(x$sign*x$vmat%*%x$y,3),
      round(x$pv,3),round(x$ci,3),round(x$tailarea,3))
    colnames(tab) = c("Var", "Stdz Coef", "P-value", "LowConfPt", "UpConfPt",
              "LowArea", "UpArea")
    rownames(tab) = rep("",nrow(tab))
    print(tab)
    
    cat(sprintf("\nEstimated stopping point from AIC rule = %i\n",x$khat))
  }

  invisible()
}

plot.fs <- function(x, breaks=TRUE, omit.zeros=TRUE) {
  
  if (x$completepath) {
    k = length(x$action)+1
    beta = cbind(x$beta,x$bls)
  } else {
    k = length(x$action)
    beta = x$beta
  }
 
  x = 1:k
  xlab = "Step"

  if (omit.zeros) {
    jj = which(rowSums(abs(beta))==0)
    if (length(jj)>0) beta = beta[jj,]
    else beta = rep(0,k)
  }

  matplot(x,t(beta),xlab=xlab,ylab="Coefficients",type="l",lty=1)
  if (breaks) abline(v=x,lty=2)
  invisible()
}


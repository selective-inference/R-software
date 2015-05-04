lar <- function(x,y, maxsteps=2000, minlam=0, verbose=FALSE,normalize=TRUE, intercept=TRUE, eps = .Machine$double.eps) {

 # We compute the least angle regression (LAR) path given
# a response vector y and predictor matrix x.  We assume
# that x has columns in general position.
    # NOTE: Rob changed X-> x, for stylistic compatibility with other functions
    
this.call=match.call()
  if (missing(y)) stop("y is missing.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 0].")
  if (missing(x)) stop("x is missing.")
  if (!is.null(x) && !is.matrix(x)) stop("x must be a matrix.")
  if (length(y)!=nrow(x)) stop("Dimensions don't match [length(y) != nrow(x)].")
  if (ncol(x) == 0) stop("There must be at least one predictor [must have ncol(x) > 0].")
  if (checkcols(x)) stop("x cannot have duplicate columns.")

  # For simplicity
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  one <- rep(1, n)
######    # added by Rob
  if (intercept) {
        meanx <- drop(one %*% x)/n
        x <- scale(x, meanx, FALSE)
        mu <- mean(y)
        y <- drop(y - mu)
    }
    else {
        meanx <- rep(0,p)
        mu <- 0
        y <- drop(y)
    }

  if (normalize) {
        normx <- sqrt(drop(one %*% (x^2)))
        nosignal <- normx/sqrt(n) < eps
        if (any(nosignal)) {
            ignores <- im[nosignal]
            inactive <- im[-ignores]
            normx[nosignal] <- eps * sqrt(n)
            if (trace) 
                cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < eps; dropped for good\n")
        }
        else ignores <- NULL
        names(normx) <- NULL
        x <- scale(x, FALSE, normx)
    }
    else {
        normx <- rep(1, p)
        ignores <- NULL
    }
#####
  # Find the first variable to enter and its sign
  uhat = t(x)%*%y
  ihit = which.max(abs(uhat))   # Hitting coordinate
  hit = abs(uhat[ihit])         # Critical lambda
  s = Sign(uhat[ihit])          # Sign
  hit = abs(uhat[ihit])

  if (verbose) {
    cat(sprintf("1. lambda=%.3f, adding variable %i, |A|=%i...",
                hit,ihit,1))
  }

  # Now iteratively find the new LAR estimate, and
  # the next critical lambda

  # Things to keep track of, and return at the end
  buf = min(maxsteps,500)
  lambda = numeric(buf)      # Critical lambdas
  action = numeric(buf)      # Action taken
  df = numeric(buf)          # Degrees of freedom
  beta = matrix(0,p,buf)     # LAR estimates
  
  lambda[1] = hit
  action[1] = ihit
  df[1] = 1
  beta[,1] = 0

  # Gamma matrix!
  Gamma = matrix(0,0,n)
  if (p>1) Gamma = rbind(Gamma,t(s*x[,ihit]+x[,-ihit]),t(s*x[,ihit]-x[,-ihit]))
  Gamma = rbind(Gamma,t(s*x[,ihit]))
  nk = nrow(Gamma)

  # M plus
  if (p>1) {
    c = t(as.numeric(Sign(t(x)%*%y)) * t(x))
    ratio = t(c[,-ihit])%*%c[,ihit]/sum(c[,ihit]^2)
    ip = 1-ratio > 0
    crit = (t(c[,-ihit])%*%y - ratio*sum(c[,ihit]*y))/(1-ratio)
    mp = max(max(crit[ip]),0)
  }
  else mp = 0
  
  # Other things to keep track of, but not return
  r = 1                      # Size of active set
  A = ihit                   # Active set
  I = Seq(1,p)[-ihit]        # Inactive set
  X1 = x[,ihit,drop=FALSE]   # Matrix X[,A]
  X2 = x[,-ihit,drop=FALSE]  # Matrix X[,I]
  k = 2                      # What step are we at?

  # Compute a skinny QR decomposition of X1
  junk = qr(X1)
  Q = qr.Q(junk,complete=TRUE)
  Q1 = Q[,1,drop=FALSE];
  Q2 = Q[,-1,drop=FALSE]
  R = qr.R(junk)
  
  # Throughout the algorithm, we will maintain
  # the decomposition X1 = Q1*R. Dimenisons:
  # X1: n x r
  # Q1: n x r
  # Q2: n x (n-r)
  # R:  r x r
    
  while (k<=maxsteps && lambda[k-1]>=minlam) {
    ##########
    # Check if we've reached the end of the buffer
    if (k > length(lambda)) {
      buf = length(lambda)
      lambda = c(lambda,numeric(buf))
      action = c(action,numeric(buf))
      df = c(df,numeric(buf))
      beta = cbind(beta,matrix(0,p,buf))
    }

    # Key quantities for the hitting times
    a = backsolve(R,t(Q1)%*%y)
    b = backsolve(R,backsolve(R,s,transpose=TRUE))
    aa = as.numeric(t(X2) %*% (y - X1 %*% a))
    bb = as.numeric(t(X2) %*% (X1 %*% b))
    
    # If the inactive set is empty, nothing will hit
    if (r==min(n,p)) hit = 0

    # Otherwise find the next hitting time
    else {
      shits = Sign(aa)
      hits = aa/(shits-bb)

      # Make sure none of the hitting times are larger
      # than the current lambda 
      hits[hits>lambda[k-1]] = 0
        
      ihit = which.max(hits)
      hit = hits[ihit]
      shit = shits[ihit]
    }

     # Stop if the next critical point is negative
    if (hit<=0) break
    
    # Record the critical lambda and solution
    lambda[k] = hit
    action[k] = I[ihit]
    df[k] = r+1
    beta[A,k] = a-hit*b
        
    # Gamma matrix!
    X2perp = X2 - X1 %*% backsolve(R,t(Q1)%*%X2)
    c = t(t(X2perp)/(shits-bb))
    Gamma = rbind(Gamma,shits*t(X2perp))
    if (ncol(c)>1) Gamma = rbind(Gamma,t(c[,ihit]-c[,-ihit]))
    Gamma = rbind(Gamma,t(c[,ihit]))
    nk = c(nk,nrow(Gamma))

    # M plus
    if (ncol(c)>1) {
      ratio = t(c[,-ihit])%*%c[,ihit]/sum(c[,ihit]^2)
      ip = 1-ratio > 0
      crit = (t(c[,-ihit])%*%y - ratio*sum(c[,ihit]*y))/(1-ratio)
      mp = c(mp,max(max(crit[ip]),0))
    }
    else mp = c(mp,0)
    
    # Update all of the variables
    r = r+1
    A = c(A,I[ihit])
    I = I[-ihit]
    s = c(s,shit)
    X1 = cbind(X1,X2[,ihit])
    X2 = X2[,-ihit,drop=FALSE]

    # Update the QR decomposition
    junk = updateQR(Q1,Q2,R,X1[,r])
    Q1 = junk$Q1
    Q2 = junk$Q2
    R = junk$R
     
    if (verbose) {
      cat(sprintf("\n%i. lambda=%.3f, adding variable %i, |A|=%i...",
                  k,hit,A[r],r))
    }
            
    # Step counter
    k = k+1
  }

  # Trim
  lambda = lambda[Seq(1,k-1)]
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
  }

  # If we reached the minimum lambda
  else if (lambda[k-1]<minlam) {
    if (verbose) {
      cat(sprintf("\nReached the minimum lambda (%.3f),",minlam))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  }
  
  # Otherwise, note that we completed the path
  else completepath = TRUE

  if (verbose) cat("\n")

beta=cbind(beta, ginv(x)%*%y) #ROB ADDED for compatibility with lars;
beta=t(beta)# for compatibility with LARS function

  out=list(lambda=lambda,action=action,df=df,beta=beta,
              completepath=completepath,Gamma=Gamma,nk=nk,mp=mp,mu=mu,meanx=meanx,normx=normx,normalize=normalize,intercept=intercept,type="LAR", call=this.call)
  class(out)="lar"
  return(out)
}


predict.lar=
function (object, newx, s, type = c("fit", "coefficients"), mode = c("step", 
    "fraction", "norm", "lambda"), ...) 
{
    mode <- match.arg(mode)
    type <- match.arg(type)

   
    if (missing(newx) & type == "fit") {
        warning("Type=fit with no newx argument; type switched to coefficients")
        type <- "coefficients"
    }
    betas <- object$beta
    
    if (object$type != "LASSO" && mode %in% c("fraction", "norm")) 
        betas = betabreaker(object)
    dimnames(betas) = list(NULL, dimnames(betas)[[2]])
    sbetas <- scale(betas, FALSE, 1/object$normx)
    kp <- dim(betas)
    k <- kp[1]
    p <- kp[2]
    steps <- seq(k)
    if (missing(s)) {
        s <- steps
        mode <- "step"
    }
    sbeta <- switch(mode, step = {
        if (any(s < 0) | any(s > k)) stop("Argument s out of range")
        steps
    }, fraction = {
        if (any(s > 1) | any(s < 0)) stop("Argument s out of range")
        nbeta <- drop(abs(sbetas) %*% rep(1, p))
        nbeta/nbeta[k]
    }, norm = {
        nbeta <- drop(abs(sbetas) %*% rep(1, p))
        if (any(s > nbeta[k]) | any(s < 0)) stop("Argument s out of range")
        nbeta
    }, lambda = {
        lambdas = object$lambda
        s[s > max(lambdas)] = max(lambdas)
        s[s < 0] = 0
        c(lambdas, 0)
    })
    sfrac <- (s - sbeta[1])/(sbeta[k] - sbeta[1])
    sbeta <- (sbeta - sbeta[1])/(sbeta[k] - sbeta[1])
    usbeta <- unique(sbeta)
    useq <- match(usbeta, sbeta)
    sbeta <- sbeta[useq]
    betas <- betas[useq, , drop = FALSE]
    coord <- approx(sbeta, seq(sbeta), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    newbetas <- ((sbeta[right] - sfrac) * betas[left, , drop = FALSE] + 
        (sfrac - sbeta[left]) * betas[right, , drop = FALSE])/(sbeta[right] - 
        sbeta[left])
    newbetas[left == right, ] <- betas[left[left == right], ]
    robject <- switch(type, coefficients = list(s = s, fraction = sfrac, 
        mode = mode, coefficients = drop(newbetas)), fit = list(s = s, 
        fraction = sfrac, mode = mode, fit = drop(scale(newx, 
            object$meanx, FALSE) %*% t(newbetas)) + object$mu))
    robject
}



larInf=function(x,y,larfit,sigma=NULL,compute.si=TRUE,alpha=.10,one.sided=TRUE,nsteps = min(nrow(x), ncol(x))){
    this.call=match.call()
    #main function for lar inference
       checkargs(x=x,y=y,larfit=larfit,sigma=sigma,alpha=alpha,nsteps=nsteps)
SMALL=1e-6
p=ncol(x)
    n=nrow(x)
nk=larfit$nk
     if(is.null(sigma) & p>=n){cat("Warning: p ge n; the value sigma=1 is used; you may want to estimate sigma using the estimateSigma function",fill=T); sigma=1}
  if(is.null(sigma) & n>p){
      sigma=sqrt(sum(lsfit(x,y)$res^2)/(n-p))
       cat("Standard deviation of noise estimated from mean squared residual",fill=T)
    }
     
vmm=vpp=pv=sigma.eta=rep(NA,nsteps)
ci=miscov=matrix(NA,nsteps,2)
for(k in 1:nsteps){
    mod=larfit$act[1:k]
 temp=(solve(t(x[,mod,drop=F])%*%x[,mod,drop=F])%*%t(x[,mod,drop=F]))
    temp=temp[nrow(temp),]
    bhat=sum(temp*y)
  eta=as.vector(temp)
    if(one.sided)eta=eta*sign(bhat)
    flip=(one.sided & sign(bhat)==-1)
    tt=sum(eta*y)
    A= -larfit$Gam[1:nk[k],]
    pp=nrow(A)
    b=rep(0,pp)
  
    vs = compute.vmvp(y, eta, A, b, pp)
     vpp = vs$vp
     vmm = vs$vm
   
    sigma.eta=sigma*sqrt(sum(eta^2))
       u=0  #null
pv[k]=1-rob.ptruncnorm(tt, vmm, vpp, u, sigma.eta)
if(!one.sided)  pv[k]=2*min(pv[k],1-pv[k])
  if(compute.si)
      {
          vs=list(vm=vmm,vp=vpp)
           junk=selection.int(y,eta,sigma, vs, alpha,flip=flip)
          ci[k,]=junk$ci
          miscov[k,]=junk$miscov
      }
 
}

pv.spacing=spacing.pval.asymp.list(y,larfit,nsteps,sigma=sigma)
   junk=covtest(larfit,x,y,sigma,nsteps)
    pv.cov=1-pexp(junk,1)

     forwardStopHat=forwardStop(pv,alpha)
    out=list(pv=pv,ci=ci,tailarea=miscov,vm=vmm,vp=vpp,sigma=sigma,alpha=alpha,act=larfit$act,pv.spacing=pv.spacing, pv.cov=pv.cov, forwardStopHat= forwardStopHat)
    out$call=this.call
   class(out)="larInf"
return(out)
}


print.larInf=function(x,digits = max(3, getOption("digits") - 3),...){
      cat("\nCall: ", deparse(x$call), "\n\n")
      cat(c("alpha=",x$alpha),fill=T)
      cat("",fill=T)
      nsteps=length(x$pv)
tab=cbind(1:nsteps,x$act[1:nsteps],round(x$pv,3),round(x$ci,3),round(x$tailarea,3),round(x$pv.spacing,3),round(x$pv.cov,3))
      dimnames(tab)=list(NULL,c("step","pred","exactPv","lowConfPt","upConfPt","lowArea","upArea","spacingPv","covtestPv"))
      print(tab)
        cat("",fill=T)
cat(c("Estimated standard deviation of noise sigma=", round(x$sigma,3)),fill=T)
         cat("",fill=T)
      cat("",fill=T)
cat(c("Estimated stopping point from forwardStop rule=", x$forwardStopHat),fill=T)
         cat("",fill=T)
  }


covtest=function(fitobj,x,y,sigma,nsteps = min(nrow(x), ncol(x))) {

    n = nrow(x)
    p = ncol(x)
    my = mean(y)
      betas <- fitobj$beta
    lambda.min.ratio = ifelse(nrow(x) < ncol(x), 0.1, 1e-04)
    jlist = unlist(fitobj$act)
        lamlist = c(fitobj$lambda, 0)
    nsteps.call = nsteps
    nsteps = length(jlist)
    nsteps = min(nsteps, which(lamlist == 0))
    nsteps = min(nsteps, nsteps.call)
    jlist = jlist[1:nsteps]
    cov0 = cov = sig = rep(NA, nsteps)
    yy = y - my
    for (j in 1:nsteps) {
            lambda = lamlist[j + 1]
                yhat = predict.lar(fitobj, x, s = lambda, type = "fit", mode = "lam")$fit   
            cov[j] = sum(yy * yhat)
            if (j == 1) {
                cov0[j] = 0
            }
            if (j > 1) {
                tt0 = which(betas[j, ] != 0)
             #     aa = update(fitobj, x = x[, tt0, drop = F])
          # aa=lar(y, x = x[, tt0, drop = F], normalize=fitobj$normalize)# use this once ryan adds the option to lar
                 aa=lar(x[, tt0, drop = F], y,normalize=fitobj$normalize,intercept=fitobj$intercept)
       #         aa$beta=cbind(aa$beta,lsfit(x[,tt0,drop=F],y)$coef[-1]) # NOTE these and next few  only needed until ryans fixes lar
    #           aa$meanx=colMeans(x[,tt0,drop=F]);aa$normx=rep(1,length(tt0))
    #            aa$mu=mean(y)
   #             aa$normalize=FALSE
                  yhat0 = predict.lar(aa, x[, tt0,drop=F], type = "fit", s = lambda, mode = "lam")$fit
                cov0[j] = sum(yy * yhat0)
            }}
        
   
  
    tt = ((cov - cov0)/sigma^2)
    
return(tt)
}

# Asymptotic spacing LAR p values list
spacing.pval.asymp.list= function(y,out,k,sigma=1) {
  if (length(out$lambda)==1) stop("The LAR path needs to be run for at least 2 steps.")
  pvals = numeric(k)
  for (i in 1:k) {
    v = out$Gamma[out$nk[i],]
    denom = sigma*sqrt(sum(v*v))
    if (i==1 && i<length(out$lambda)) {
      pvals[1] = (1-pnorm(out$lambda[i]/denom))/
        (1-pnorm(out$lambda[i+1]/denom))
    }
    else if (i<length(out$lambda)) {
    #  pvals[i] = (pnorm(out$lambda[i-1]/denom)-pnorm(out$lambda[i]/denom))/
     #   (pnorm(out$lambda[i-1]/denom)-pnorm(out$lambda[i+1]/denom))
        pvals[i]=1-ptruncnorm(out$lambda[i],out$lambda[i+1],out$lambda[i-1],0,denom)
    }
    else {
     # pvals[i] = (pnorm(out$lambda[i-1]/denom)-pnorm(out$lambda[i]/denom))/
      #  (pnorm(out$lambda[i-1]/denom)-pnorm(0))
         pvals[i]=1-ptruncnorm(out$lambda[i],out$lambda[i+1],0,0,denom)
    }
  }
  return(pvals)
}
 plot.lar=
function (x, xvar = c("norm", "df", "arc.length", "step"), breaks = TRUE, 
    plottype = c("coefficients", "Cp"), omit.zeros = TRUE, eps = 1e-10, 
    ...) 
{
    object <- x
    plottype <- match.arg(plottype)
    xvar <- match.arg(xvar)
    coef1 <- object$beta
    if (x$type != "LASSO" && xvar == "norm") 
        coef1 = betabreaker(x)
    stepid = trunc(as.numeric(dimnames(coef1)[[1]]))
    coef1 <- scale(coef1, FALSE, 1/object$normx)
    if (omit.zeros) {
        c1 <- drop(rep(1, nrow(coef1)) %*% abs(coef1))
        nonzeros <- c1 > eps
        cnums <- seq(nonzeros)[nonzeros]
        coef1 <- coef1[, nonzeros, drop = FALSE]
    }
    else cnums <- seq(ncol(coef1))
    s1 <- switch(xvar, norm = {
        s1 <- apply(abs(coef1), 1, sum)
        s1/max(s1)
    }, df = object$df, arc.length = cumsum(c(0, object$arc.length)), 
        step = seq(nrow(coef1)) - 1)
    xname <- switch(xvar, norm = "|beta|/max|beta|", df = "Df", 
        arc.length = "Arc Length", step = "Step")
    if (plottype == "Cp") {
        Cp <- object$Cp
        plot(s1, Cp, type = "b", xlab = xname, main = object$type, 
            ...)
    }
    else {
        matplot(s1, coef1, xlab = xname, ..., type = "b", pch = "*", 
            ylab = "Standardized Coefficients")
        title(object$type, line = 2.5)
        abline(h = 0, lty = 3)
        axis(4, at = coef1[nrow(coef1), ], labels = paste(cnums), 
            cex = 0.8, adj = 0)
        if (breaks) {
            axis(3, at = s1, labels = paste(stepid), cex = 0.8)
            abline(v = s1)
        }
    }
    invisible()
}


# Make sure that no two columms of A are the same
# (this works with probability one).

checkcols <- function(A) {
  b = rnorm(nrow(A))
  a = sort(t(A)%*%b)
  if (any(diff(a)==0)) return(TRUE)
  return(FALSE)
}

# Downdate the QR factorization, after a column has
# been deleted. Here Q1 is m x n, Q2 is m x k, and
# R is n x n.

downdateQR <- function(Q1,Q2,R,col) {
  m = nrow(Q1)
  n = ncol(Q1)
  
  a = .C("downdate1",
    Q1=as.double(Q1),
    R=as.double(R),
    col=as.integer(col-1),
    m=as.integer(m),
    n=as.integer(n),
    dup=FALSE,
    package="genlasso")

  Q1 = matrix(a$Q1,nrow=m)
  R = matrix(a$R,nrow=n)

  # Re-structure: add a column to Q2, delete one from
  # Q1, and trim R
  Q2 = cbind(Q2,Q1[,n])
  Q1 = Q1[,-n,drop=FALSE]
  R = R[-n,-col,drop=FALSE]

  return(list(Q1=Q1,Q2=Q2,R=R))
}

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
    package="polypath")

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


betabreaker=
function (object) 
{
    frac.step = function(x) {
        x[1] - 1 - x[2]/(x[3] - x[2])
    }
    beta = object$beta
    sbeta = sign(beta)
    kp = dim(beta)
    k = kp[1]
    p = kp[2]
    dsbeta = abs(sbeta[-1, ] - sbeta[-k, ])
    if (any(dsbeta == 2)) {
        bbeta = matrix(cbind(step.end = rep(1:(k - 1), p), beta.start = as.vector(beta[-k, 
            ]), beta.end = as.vector(beta[-1, ]))[dsbeta == 2], 
            ncol = 3)
        fsteps = apply(bbeta, 1, frac.step)
        new.beta = predict(object, type = "coefficient", s = fsteps + 
            1, mode = "step")$coef
        new.beta = rbind(beta, new.beta)
        fo = c(seq(k) - 1, fsteps)
        beta = new.beta[order(fo), ]
        dimnames(beta)[[1]] = format(round(sort(fo), 2))
    }
    beta
}

# We compute the least angle regression (LAR) path given
# a response vector y and predictor matrix X.  We assume
# that X has columns in general position.

lar <- function(X,y, maxsteps=2000, minlam=0, verbose=FALSE,normalize=TRUE, intercept=TRUE, eps = .Machine$double.eps) {
this.call=match.call()
  if (missing(y)) stop("y is missing.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 0].")
  if (missing(X)) stop("X is missing.")
  if (!is.null(X) && !is.matrix(X)) stop("X must be a matrix.")
  if (length(y)!=nrow(X)) stop("Dimensions don't match [length(y) != nrow(X)].")
  if (ncol(X) == 0) stop("There must be at least one predictor [must have ncol(X) > 0].")
  if (checkcols(X)) stop("X cannot have duplicate columns.")

  # For simplicity
  y = as.numeric(y)
  n = nrow(X)
  p = ncol(X)
  one <- rep(1, n)
######    # added by Rob
  if (intercept) {
        meanx <- drop(one %*% X)/n
        X <- scale(X, meanx, FALSE)
        mu <- mean(y)
        y <- drop(y - mu)
    }
    else {
        meanx <- rep(0,p)
        mu <- 0
        y <- drop(y)
    }

  if (normalize) {
        normx <- sqrt(drop(one %*% (X^2)))
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
        X <- scale(X, FALSE, normx)
    }
    else {
        normx <- rep(1, p)
        ignores <- NULL
    }
#####
  # Find the first variable to enter and its sign
  uhat = t(X)%*%y
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
  if (p>1) Gamma = rbind(Gamma,t(s*X[,ihit]+X[,-ihit]),t(s*X[,ihit]-X[,-ihit]))
  Gamma = rbind(Gamma,t(s*X[,ihit]))
  nk = nrow(Gamma)

  # M plus
  if (p>1) {
    c = t(as.numeric(Sign(t(X)%*%y)) * t(X))
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
  X1 = X[,ihit,drop=FALSE]   # Matrix X[,A]
  X2 = X[,-ihit,drop=FALSE]  # Matrix X[,I]
  k = 2                      # What step are we at?

  # Compute a skinny QR decomposition of X1
  x = qr(X1)
  Q = qr.Q(x,complete=TRUE)
  Q1 = Q[,1,drop=FALSE];
  Q2 = Q[,-1,drop=FALSE]
  R = qr.R(x)
  
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
    x = updateQR(Q1,Q2,R,X1[,r])
    Q1 = x$Q1
    Q2 = x$Q2
    R = x$R
     
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

beta=cbind(beta,lsfit(X,y)$coef[-1]) #ROB ADDED for compatibility with lars; NOTE this needs to be fixed, eg for p>n

beta=t(beta)# for compatibility with LARS function

  out=list(lambda=lambda,action=action,df=df,beta=beta,
              completepath=completepath,Gamma=Gamma,nk=nk,mp=mp,mu=mu,meanx=meanx,normx=normx,normalize=normalize,intercept=intercept,type="LAR", call=this.call)
  class(out)="lars"
  return(out)
}

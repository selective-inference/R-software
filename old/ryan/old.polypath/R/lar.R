# We compute the least angle regression (LAR) path given
# a response vector y and predictor matrix X.  We assume
# that X has columns in general position.

lar <- function(y, X, maxsteps=2000, minlam=0, verbose=FALSE) {
 
  if (missing(y)) stop("y is missing.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 1].")
  if (missing(X)) stop("X is missing.")
  if (!is.null(X) && !is.matrix(X)) stop("X must be a matrix.")
  if (length(y)!=nrow(X)) stop("Dimensions don't match [length(y) != nrow(X)].")
  if (checkcols(X)) stop("X cannot have duplicate columns.")

  # For simplicity
  y = as.numeric(y)
  n = nrow(X)
  p = ncol(X)

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
  Gamma = rbind(
    t(s*X[,ihit]+X[,Seq(1,p)[-ihit]]),
    t(s*X[,ihit]-X[,Seq(1,p)[-ihit]]),
    t(s*X[,ihit]))
  nk = nrow(Gamma)
  
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
    c = t(rbind(t(X2perp)/(1-bb),t(X2perp)/(-1-bb)))
    c0 = Gamma[nrow(Gamma),]
    imax = ihit + (p-r)*(shit==-1)
    ii = t(c) %*% y <= lambda[k-1]
    jj = ii & 1:(2*(p-r))!=imax

    if (any(ii)) Gamma = rbind(Gamma,t(c0-c[,ii]))
    if (any(!ii)) Gamma = rbind(Gamma,t(c[,!ii]-c0))
    if (any(jj)) Gamma = rbind(Gamma,t(c[,imax]-c[,jj]))
    Gamma = rbind(Gamma,t(c[,imax]))
    nk = c(nk,nrow(Gamma))
    
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

  return(list(lambda=lambda,action=action,df=df,beta=beta,
              completepath=completepath,Gamma=Gamma,nk=nk))
}

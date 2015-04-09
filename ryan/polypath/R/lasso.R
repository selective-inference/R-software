# We compute the lasso solution path given a response vector y
# and predictor matrix X.  This is:
# \min ||y - X\beta||_2^2 + \lambda ||\beta||_1
# We assume that X has columns in general position.

lasso <- function(y, X, maxsteps=2000, minlam=0, verbose=FALSE) {
 
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

  # Gamma matrix!
  Gamma = matrix(0,0,n)
  if (p>1) Gamma = rbind(Gamma,t(s*X[,ihit]+X[,-ihit]),t(s*X[,ihit]-X[,-ihit]))
  Gamma = rbind(Gamma,t(s*X[,ihit]))
  nk = nrow(Gamma)
  #nkk = nk
  #sb = hit
  
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
  # the decomposition X1 = Q1*R. Dimensions:
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

    # Key quantities for hitting and leaving times
    a = backsolve(R,t(Q1)%*%y)
    b = backsolve(R,backsolve(R,s,transpose=TRUE))
    aa = as.numeric(t(X2) %*% (y - X1 %*% a))
    bb = as.numeric(t(X2) %*% (X1 %*% b))
    
    # If the inactive set is empty, nothing will hit
    ### NOTE: checking r==min(n,p), because of the Q2 factor
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

    # If nothing is on the boundary, nothing will leave
    if (r==0) leave = 0
      
    # Otherwise, find the next leaving time
    else {
      leaves = a/b
        
      # Make sure none of the leaving times are larger
      # than the current lambda 
      leaves[leaves>lambda[k-1]] = 0

      # If a variable just entered, then make sure it 
      # cannot leave
      if (action[k-1]>0) leaves[r] = 0
      
      ileave = which.max(leaves)
      leave = leaves[ileave]
    }

    ##########
    # Stop if the next critical point is negative
    if (hit<=0 && leave<=0) break
    
    # If a hitting time comes next
    if (hit > leave) {
      # Record the critical lambda and solution
      lambda[k] = hit
      action[k] = I[ihit]
      df[k] = r+1
      beta[A,k] = a-hit*b
      
      # Gamma matrix!
      c0 = Gamma[nrow(Gamma),]
      # Hitting times
      X2perp = X2 - X1 %*% backsolve(R,t(Q1)%*%X2)
      chit = t(t(X2perp)/(shits-bb))
      Gamma = rbind(Gamma,shits*t(X2perp))
      if (ncol(chit)>1) Gamma = rbind(Gamma,t(chit[,ihit]-chit[,-ihit]))
      # Leaving times
      X1pinv = backsolve(R,t(Q1))
      clea = t(X1pinv/b)
      if (action[k-1]>0) clea = clea[,-r,drop=FALSE]
      ii.lea = (t(clea) %*% y <= lambda[k-1])
      if (any(ii.lea)) Gamma = rbind(Gamma,t(c0-clea[,ii.lea]))
      if (any(!ii.lea)) Gamma = rbind(Gamma,t(clea[,!ii.lea]-c0))
      ### NOTE: it's important to redefine both of clea.s and ilea.s,
      ### since we only want to look at variables whose leaving time
      ### is <= lambda[k-1]
      clea.s = clea[,ii.lea,drop=FALSE]
      ilea.s = which.max(t(clea.s)%*%y)
      if (ncol(clea.s)>1) Gamma = rbind(Gamma,t(clea.s[,ilea.s]-clea.s[,-ilea.s]))
      # Hitting time comes next
      if (ncol(clea.s)>0) Gamma = rbind(Gamma,t(chit[,ihit]-clea.s[,ilea.s]))
      Gamma = rbind(Gamma,t(chit[,ihit]))
      nk = c(nk,nrow(Gamma))
      #nkk = c(nkk,nrow(Gamma)-1)
      #sb = c(sb,t(cmax.hit)%*%y)
      
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
    }

    # Otherwise a leaving time comes next
    else {
      # Record the critical lambda and solution
      lambda[k] = leave
      action[k] = -A[ileave]
      df[k] = r-1
      beta[A,k] = a-leave*b
         
      # Gamma matrix!
      c0 = Gamma[nrow(Gamma),]
      # Hitting times
      X2perp = X2 - X1 %*% backsolve(R,t(Q1)%*%X2)
      if (r==min(n,p)) chit = matrix(0,n,0)
      else chit = t(t(X2perp)/(shits-bb))
      if (ncol(chit)>0) Gamma = rbind(Gamma,shits*t(X2perp))
      if (ncol(chit)>1) Gamma = rbind(Gamma,t(chit[,ihit]-chit[,-ihit]))
      # Leaving times
      X1pinv = backsolve(R,t(Q1))
      clea = t(X1pinv/b)
      if (action[k-1]>0) clea = clea[,-r,drop=FALSE]
      ii.lea = (t(clea) %*% y <= lambda[k-1])
      Gamma = rbind(Gamma,t(c0-clea[,ii.lea]))
      if (any(!ii.lea)) Gamma = rbind(Gamma,t(clea[,!ii.lea]-c0))
      ### NOTE: it's important to redefine both of clea.s and ilea.s,
      ### since we only want to look at variables whose leaving time
      ### is <= lambda[k-1]
      clea.s = clea[,ii.lea,drop=FALSE]
      ilea.s = which.max(t(clea.s)%*%y)
      Gamma = rbind(Gamma,t(clea.s[,ilea.s]-clea.s[,-ilea.s]))
      # Leaving time comes next
      if (ncol(chit)>0) Gamma = rbind(Gamma,t(clea.s[,ilea.s]-chit[,ihit]))
      Gamma = rbind(Gamma,t(clea.s[,ilea.s]))
      nk = c(nk,nrow(Gamma))
      #nkk = c(nkk,nrow(Gamma)-1)
      #sb = c(sb,t(cmax.hit)%*%y)
     
      # Update all of the variables
      r = r-1
      I = c(I,A[ileave])
      A = A[-ileave]
      s = s[-ileave]
      X2 = cbind(X2,X1[,ileave])
      X1 = X1[,-ileave,drop=FALSE]

      # Downdate the QR decomposition
      x = downdateQR(Q1,Q2,R,ileave)
      Q1 = x$Q1
      Q2 = x$Q2
      R = x$R

      if (verbose) {
        cat(sprintf("\n%i. lambda=%.3f, deleting variable %i, |A|=%i...",
                    k,leave,I[p-r],r))
      }
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

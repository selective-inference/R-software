## Psub <- x %*% ginv(x)
## Pfull <- X %*% ginv(X)
## Pz <- diag(rep(1,n)) - Psub
## PM <- diag(rep(1,n)) - Pfull
## R1 <- Pz %*% y
## R2 <- PM %*% y
## z <- y - R1
## C <- sum(diag(Pfull-Psub))/sum(diag(Psub))
## norm2R1 <- sum(R1^2)
## norm2R2 <- sum(R2^2)
## TF <- (norm2R1-norm2R2)/(C*norm2R2)
## r = sqrt(norm2R2)
## Vdelta <- (R1-R2)/sqrt(norm2R1-norm2R2)
## V2 <- R2/r

# -----------------------------------------------------------

#' Form univariate selection interval for a given contrast
#'
#' @param obj Object returned by \link{groupfs} function
#' @param x Design matrix
#' @param y Outcome
#' @param index index of variable grouping
#' @return An interval or union of intervals
groupfsInf <- function(obj, sigma = NULL, normalize = TRUE, projs = NULL, tol = 1e-15) {

  x <- obj$x
  n <- nrow(x)
  p <- ncol(x)
  y <- obj$y
  active <- obj$action
  maxsteps <- attr(obj, "maxsteps")
  sp <- attr(obj, "sp")
  index <- obj$index
  Eindex <- which(index %in% active)
  Ep <- length(Eindex)
  
  pvals <- numeric(maxsteps)

  if (is.null(sigma)) {
    if (n >= 2*p) {
      sigma <- sqrt(sum(lsfit(x, y, intercept = obj$intercept)$res^2)/(n-p-obj$intercept))
    } else {
      #sigma = sd(y)
      sigma = sqrt(obj$log$RSS[length(obj$log$RSS)]/(n-Ep-obj$intercept))
      warning(paste(sprintf("p > n/2, and sigmahat = %0.3f used as an estimate of sigma;",sigma), "you may want to use the estimateSigma function"))
    }
  }
  
  # Compute p-value for each active group
  for (j in 1:maxsteps) {
    i <- active[j]
    #Ug <- fit$maxprojs[[j]]

    # Form projection onto active set \i
    # and project x_i orthogonally
    x_i <- x[,which(index == i), drop = FALSE]      
    if (length(active) > 1) {
        minus_i <- setdiff(active, i)
        x_minus_i <- svd(x[,which(index %in% minus_i), drop = FALSE])$u
        x_i <- x_i - x_minus_i %*% t(x_minus_i) %*% x_i    
    }

    # Project y onto what remains of x_i
    Ugtilde <- svd(x_i)$u
    R <- t(Ugtilde) %*% y
    # R ~ N(0, \sigma^2 P) under H_0 : P\mu = 0
    # Z is independent of R, we effectively condition on Z
    # Unit vector used for finding the truncation interval
          # Fix scale
          # What is conditional variance?
    scale <- sigma
    df <- ncol(Ugtilde)
    # The truncated chi
    TC <- sqrt(sum(R^2))
    
    # For each step...
    L <- check_inequalities(obj, x, y, index, sp, TC, R, Ugtilde)

    # Any additional projections, e.g. from cross-validation?
    if (!is.null(projs)) L <- c(L, projs)
    
    # Compute intersection:
    E <- interval_complement(do.call(interval_union, L), check_valid = FALSE)
    
    # E is now potentially a union of intervals
    if (length(E) == 0) {
#        print(L)
        stop("Trivial intersection")
    }

    # Sum truncated cdf over each part of E
    denom <- do.call(sum, lapply(1:nrow(E), function(v) {
      tchi_interval(E[v,1], E[v,2], scale, df)
    }))
    
    # Sum truncated cdf from observed value to max of
    # truncation region
    numer <- do.call(sum, lapply(1:nrow(E), function(v) {
      lower <- E[v,1]
      upper <- E[v,2]
      if (upper > TC) {
        # Observed value is left of this interval's right endpoint
          
        if (lower < TC) {
          # Observed value is in this interval
          return(tchi_interval(TC, upper, scale, df))
          
        } else {
          # Observed value is not in this interval
          return(tchi_interval(lower, upper, scale, df))
        }
        
      } else {
        # Observed value is right of this entire interval
        return(0)
      }
    }))

    # Survival function
    value <- numer/denom
    if (is.nan(value)) {
      value <- 0
      warning("P-value NaN of the form 0/0 converted to 0.")
    }
    # Force p-value to lie in the [0,1] interval
    # in case of numerical issues 
    # value <- max(0, min(1, value))
    pvals[j] <- value
  }
  names(pvals) <- obj$action
  invisible(pvals)
}


# -----------------------------------------------------------

tchi_interval <- function(lower, upper, scale, df) {
  a <- (lower/scale)^2
  b <- (upper/scale)^2    
  if (b == Inf) {
      integral <- pchisq(a, df, lower.tail = FALSE)
  } else {
      integral <- pchisq(b, df) - pchisq(a, df)
  }
  if ((integral < .Machine$double.eps) && (b < Inf)) {
      integral <- num_int_chi(a, b, df)
  }
  return(integral)
}

num_int_chi <- function(a, b, df, nsamp = 10000) {
  grid <- seq(from=a, to=b, length.out=nsamp)
  integrand <- dchisq(grid, df)
  return((b-a)*mean(integrand))
}


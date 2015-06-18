# Forward stepwise implementation with groups and truncated \chi p-values

# -----------------------------------------------------------

#' Forward stepwise
#'
#' @param x Design matrix
#' @param y Response vector
#' @param index Group membership indicator of length p
#' @param Sigma Error covariance matrix
#' @param steps Maximum number of steps for forward stepwise
#' @param normalize Should the design matrix be normalized?
#' @return Index \code{imax} of added group, residualized data, truncated chi p-value
fstep <- function(x, y, index, Sigma, steps, normalize = TRUE, k = 0, ...) UseMethod("fstep")

fstep.default <- function(x, y, index, Sigma, steps, normalize = TRUE, k = 0, verbose=FALSE, ...) {

  p <- ncol(x)
  n <- nrow(x)

  # Assume no groups if index not specified
  if (missing(index)) index <- 1:p

  # Group labels
  labels <- unique(index)
  inactive <- labels
  active <- c()
  
  if (missing(Sigma)) {
    Sigma <- 1
    warning("Error covariance unspecified, using diag(1)")
  }
  
  if (missing(steps)) steps <- min(n, p) - 1
  
  # Form whitening matrix Sigma^(-1/2)
  if (is.matrix(Sigma)) {
    svdSigma <- svd(Sigma)
    whitener <- svdSigma$u %*% diag(1/sqrt(svdSigma$d)) %*% t(svdSigma$u)
  }
  
  # Initialize copies of data for loop
  y.update <- y - mean(y)
  y.last <- y.update  
  x.update <- x

  # Center and scale design matrix
  if (normalize) {
    x.update <- scale_groups(x.update, index)
  }

  # Store all projections computed along the path
  projections <- list()

  # Store other information from each step
  path.info <- data.frame(imax=integer(), k=integer(), L=numeric(), RSS=numeric(), RSSdrop=numeric(), chisq=numeric(), projections=numeric())

  # Begin main loop
  for (step in 1:steps) {
      
    added <- add1.fstep(x.update, y.update, index, labels, inactive, k, Sigma)
    
    # Group to be added
    imax <- added$imax

    # Projection to added group
    P.imax <- added$projections[[as.character(imax)]]
    
    inactive <- setdiff(inactive, imax)
    active <- union(active, imax)    
    inactive.inds <- which(!index %in% active)

    # Rank of group
    added$k <- sum(diag(P.imax))

    # Regress added group out of y and inactive x
    y.update <- y.update - P.imax %*% y.update
    x.update[, inactive.inds] <- x.update[, inactive.inds] - P.imax %*% x.update[, inactive.inds]

    # Form projections P_{imax} - P_{inactive}
    # These are the quadratic forms Q characterizing imax
    # by t(y) %*% Q %*% y >= 0
    projs <- lapply(names(added$projections), function(i) {
      return(P.imax - added$projections[[i]])
    })
    names(projs) <- names(added$projections)
    # Delete the 0 element formed by P_imax - P_imax
    # and store the rest
    projections[[step]] <- projs[-which(names(projs) == as.character(imax))]

    # Compute RSS for unadjusted chisq p-values
    if (is.matrix(Sigma)) {
      added$RSS <- t(y.update) %*% whitener %*% y.update
      scale.chisq <- 1
    } else {
      added$RSS <- sum(y.update^2)
      scale.chisq <- Sigma
    }
    added$RSSdrop <- sum((y.last - y.update)^2)
    added$chisq <- pchisq(added$RSSdrop/scale.chisq, lower.tail=FALSE, df = added$k)
    y.last <- y.update

    # Projections are stored separately
    added$projections <- NULL
    path.info <- rbind(path.info, data.frame(added))
    
    if (verbose) print(added)
  }

  # Create output object
  value <- list(variable=path.info$imax, L=path.info$L, projections = projections, log = path.info)
  class(value) <- "fstep"
  attr(value, "n") <- nrow(x)
  attr(value, "p") <- ncol(x)
  attr(value, "labels") <- labels
  attr(value, "index") <- index
  attr(value, "steps") <- steps
  if (!is.null(attr(x, "varnames"))) {
    attr(value, "varnames") <- attr(x, "varnames")
  } else {
    attr(value, "varnames") <- colnames(x)
  }
  if (!is.matrix(Sigma)) attr(value, "scale") <- Sigma

  invisible(value)
}

# -----------------------------------------------------------

#' Add one group to the model
#'
#' @param x Design matrix
#' @param y Response vector
#' @param index Group membership indicator of length p
#' @param inactive Indices of inactive groups
#' @param Sigma Error covariance matrix
#' @return Index \code{imax} of added group, residualized data, truncated chi p-value
add1.fstep <- function(x, y, index, labels, inactive, k, Sigma, ...) {

  n <- nrow(x)
  p <- ncol(x)

  if (missing(Sigma)) {
    Sigma <- 1
    warning("Error covariance unspecified, using diag(1)")
  }

  # Use characters to avoid issues where
  # list() populates NULL lists in the positions
  # the active variables
  ### Question for later: does this slow down lapply?
  keys = as.character(inactive)  

  # Compute projections: X_i %*% X_i^\dagger
  projections <- lapply(keys, function(i) {
    inds <- which(index == i)
    xi <- x[,inds]
    return(xi %*% ginv(xi))
  })

  names(projections) <- keys

  # Compute sums of squares to determine which group is added
  # penalized by rank of group if k > 0
  terms <- lapply(keys, function(i) {
    t(y) %*% projections[[i]] %*% y - ifelse(k > 0, k * sum(diag(projections[[i]])), 0)
  })

  # Maximizer = group to be added
  terms.maxind <- which.max(terms)
  imax <- inactive[terms.maxind]
  maxinds <- which(index == imax)
  L <- terms[[terms.maxind]]

  return(list(imax=imax, L=L, projections = projections))
}


# -----------------------------------------------------------

#' Form univariate selection interval for a given contrast
#'
#' @param fit fitted model, the result of \link{fstep}
#' @param x Design matrix
#' @param y Outcome
#' @param index index of variable grouping
#' @return An interval or union of intervals
interval.fstep <- function(fit, x, y, index, normalize = TRUE) {

  projs <- fit$projections
  pvals <- c()
  
  if (normalize) {
    x <- scale_groups(x, index)
  }  

  # Compute p-value for each active group
  for (i in fit$variable) {

    # Form projection onto active set \i
    # and project x_i orthogonally
    minus_i <- setdiff(fit$variable, i)
    x_minus_i <- x[,which(index %in% minus_i), drop = FALSE]
    x_i <- x[,which(index == i), drop = FALSE]
    x_i <- x_i - x_minus_i %*% ginv(x_minus_i) %*% x_i

    # Project y onto what remains of x_i
    P <- x_i %*% ginv(x_i)
    # R ~ N(0, \sigma^2 P) under H_0 : P\mu = 0
    R <- P %*% y
    # Z is independent of R, we effectively condition on Z
    Z <- y - R
    # Unit vector used for finding the truncation interval
    U <- R / sqrt(sum(R^2))
          # Fix scale
          # What is conditional variance?
    scale <- 1
    k <- sum(diag(P))
    # The truncated chi
    TC <- sqrt(sum(R^2))
    
    # For each step...
    L <- lapply(1:length(projs), function(j) {
        Q <- projs[[j]]

        # For each inactive group at that step...
        LL <- lapply(names(Q), function(l) {

            # FIX: which terms are shifted when k is nonzero?
            # ifelse(k > 0, k * sum(diag(projections[[i]])), 0)

        # The quadratic form corresponding to
        # (t*U + Z)^T %*% Q[[l]] %*% (t*U + Z) \geq 0
        # we find the roots in t, if there are any
        # and return the interval of potential t
        a = t(U) %*% Q[[l]] %*% U
        b = 2 * t(U) %*% Q[[l]] %*% Z
        c = t(Z) %*% Q[[l]] %*% Z
        disc <- b^2 - 4*a*c
        b2a <- -b/(2*a)

        if (disc > 0) {
          # Real roots
          pm <- sqrt(disc)/(2*a)
          endpoints <- c(b2a - pm, b2a + pm)
          
        } else {
            
          # No real roots
          if (a > 0) {
            # Quadratic form always positive
            return(Intervals(c(0,Inf)))
          } else {
            # Quadratic form always negative
            stop("Infeasible!")
          }
        }
          
        if (a > 0) {
          # Parabola opens upward
            
          if (min(endpoints) > 0) {
            # Both roots positive, union of intervals
            return(Intervals(rbind(c(0, min(endpoints)), c(max(endpoints), Inf))))
            
          } else {
            # At least one negative root
            return(Intervals(c(max(0, max(endpoints)), Inf)))
          }
          
        } else {
          # Parabola opens downward
          if (max(endpoints) < 0) {
              
            # Positive quadratic form only when t negative
            stop("Error: infeasible")
          } else {
              
            # Part which is positive
            return(Intervals(c(max(0, min(endpoints)), max(endpoints))))
          }
        }
      })
     # LL is a list of intervals
     return(LL)
    })
    
    # L is now a list of lists of intervals
    # Compute intersection:
    E <- do.call(interval_intersection, unlist(L, recursive = F))
    
    # E is now potentially a union of intervals
    if (length(E) == 0) stop("Trivial intersection")

    # Sum truncated cdf over each part of E
    denom <- do.call(sum, lapply(1:nrow(E), function(v) {
      tchi_interval(E[v,1], E[v,2], scale, k)
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
          return(tchi_interval(TC, upper, scale, k))
          
        } else {
          # Observed value is not in this interval
          return(tchi_interval(lower, upper, scale, k))
        }
        
      } else {
        # Observed value is right of this entire interval
        return(0)
      }
    }))

    # Survival function
    pvals <- c(pvals, numer/denom)    
  }
  names(pvals) <- fit$variable
  invisible(pvals)
}

# -----------------------------------------------------------

tchi_interval <- function(lower, upper, scale, df) {
    if (upper == Inf) {
      return(pchisq(lower^2/scale^2, df, lower.tail = FALSE))
    } else {
      return(pchisq(upper^2/scale^2, df) - pchisq(lower^2/scale^2, df))
    }
}


#' Center and scale design matrix by groups
#'
#' @param x Design matrix
#' @param index Group membership indicator of length p
#' @param center Center groups, default TRUE
#' @param scale Scale groups by Frobeniusm norm, default TRUE
#' @return Scaled design matrix
scale_groups <- function(x, index, center = TRUE, scale = TRUE, ...) UseMethod("scale_groups")

scale_groups.default <- function(x, index, center = TRUE, scale = TRUE, ...) {

  for (j in unique(index)) {
    inds <- which(index == j)
    if (center) x[, inds] <- x[, inds] - mean(x[, inds])
    if (scale) {
      normsq <- sum(x[, inds]^2)
      if (normsq > 0) x[, inds] <- x[, inds] / sqrt(normsq)
    }

  }
  return(x)
}


# -----------------------------------------------------------

# -----------------------------------------------------------

# ----------------- Ignore code below here ------------------

# -----------------------------------------------------------

# -----------------------------------------------------------




#' Compute AIC and return active set of model minimizing the criterion
#'
#' @param fit fitted model, the result of \link{fstep}
#' @param scale Known/estimated scaling parameter
#' @param k Multiplier of model size, use \code{k = 2} for AIC, \code{k = log(n)} for BIC, or \code{k = 2log(p)} for RIC.
#' @return Matrix with two columns, edf and aic of fit, with rows indexed by step
extractAIC.fstep <- function(fit, scale, k = 2, ...) {

  if (class(fit) != "fstep") stop("Incorrect object type of fit")

  if (missing(scale)) {
    scale <- attributes(fit)$scale
    if (is.null(scale)) {
      warning("Missing scale, assuming scale = 1")
      scale <- 1
    }
  }

  index <- attr(fit, "index")
  edf <- cumsum(fit$log$k)
  aic <- fit$log$RSS + k * edf

  return(cbind(edf, aic))
}

# -----------------------------------------------------------

#' Determine active sets using various stopping rules
#'
#' @param fit fitted model, the result of \link{fstep}
#' @return Matrix with two columns, edf and aic of fit, with rows indexed by step
model_select <- function(fit, alpha = .1, ...) {

  if (class(fit) != "fstep") stop("Incorrect object type of fit")

  # Avoid integer(0) issues with which()
  max_which <- function(l) {
    wl <- which(l)
    if (any(wl)) return(max(wl))
    return(NULL)
  }

  models <- list()

  n <- attr(fit, "n")
  p <- attr(fit, "p")

  # Last rule
  ind <- max_which(fit$p.value < alpha)
  if (is.null(ind)) {
    models$last <- 0
  } else {
    models$last <- fit$variable[1:ind]
  }

  # ForwardStop
  renyi.transform <- - cumsum(log(1-fit$p.value)) / 1:length(fit$p.value)
  ind <- max_which(renyi.transform <= alpha)
  if (is.null(ind)) {
    models$forward <- 0
  } else {
    models$forward <- fit$variable[1:ind]
  }

  # RIC
  ind <- which.min(extractAIC(fit, k = 2*log(p))[, 2])
  models$RIC <- fit$variable[1:ind]

  # BIC
  ind <- which.min(extractAIC(fit, k = log(n))[, 2])
  models$BIC <- fit$variable[1:ind]

  return(models)
}

# -----------------------------------------------------------

#' Plot \code{T\chi} and \code{\chi^2} p-values of a fit
#'
#' @param fit fitted model, the result of \link{fstep}
plot.fstep <- function(fit, ...) {
  xrange <- 1:attr(fit, "steps")
  plot.default(x = xrange,  y = fit$p.value, xlab = "Variable", ylab = "P-value", cex = .5, xaxt = "n", pch = 19, ylim = c(0,1))
  points.default(x = xrange, y = fit$log$chisq, cex = .5)
  varnames <- attr(fit, "varnames")[fit$variable]
  if (!is.null(varnames)) {
    axis(1, at = xrange, labels=FALSE)
    text(x = xrange, y=-.08, labels = varnames, srt = 45, pos = 1, xpd = TRUE, cex = .4)
  }
}


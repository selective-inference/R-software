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
groupfs <- function(x, y, index, Sigma, steps, normalize = TRUE, k = 0, ...) UseMethod("groupfs")

groupfs.default <- function(x, y, index, Sigma, steps, normalize = TRUE, k = 0, verbose=FALSE, ...) {

  p <- ncol(x)
  n <- nrow(x)

  if (steps >= n) steps <- n-1

  # Assume no groups if index not specified
  if (missing(index)) index <- 1:p

  # Group labels
  labels <- unique(index)
  inactive <- labels
  active <- c()
  
  if (missing(Sigma)) {
    Sigma <- 1
    #warning("Error covariance unspecified, using diag(1)")
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

  modelrank <- 1
  
  # Begin main loop
  for (step in 1:steps) {
      
    added <- add1.groupfs(x.update, y.update, index, labels, inactive, k, Sigma)
    
    # Group to be added
    imax <- added$imax

    # Projection to added group
    P.imax <- added$projections[[as.character(imax)]]
    
    inactive <- setdiff(inactive, imax)
    active <- union(active, imax)    
    inactive.inds <- which(!index %in% active)

    # Rank of group
    added$k <- sum(diag(P.imax))
    modelrank <- modelrank + added$k

    # Stop without adding if model has become saturated
    if (modelrank >= n) break

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
  class(value) <- "groupfs"
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
add1.groupfs <- function(x, y, index, labels, inactive, k, Sigma, ...) {

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
#' @param fit fitted model, the result of \link{groupfs}
#' @param x Design matrix
#' @param y Outcome
#' @param index index of variable grouping
#' @return An interval or union of intervals
interval.groupfs <- function(fit, x, y, index, k = 0, normalize = TRUE, tol = 1e-15) {

  pvals <- c()
  y <- y - mean(y)
  
  if (normalize) {
    x <- scale_groups(x, index)
  }  

  active <- fit$variable
  
  # Compute p-value for each active group
  for (i in active) {

    # Form projection onto active set \i
    # and project x_i orthogonally
    x_i <- x[,which(index == i), drop = FALSE]      
    if (length(active) > 1) {
        minus_i <- setdiff(active, i)
        x_minus_i <- x[,which(index %in% minus_i), drop = FALSE]
        x_i <- x_i - x_minus_i %*% ginv(x_minus_i) %*% x_i    
    }

    # Project y onto what remains of x_i
    P <- x_i %*% ginv(x_i)
    R <- P %*% y
    # R ~ N(0, \sigma^2 P) under H_0 : P\mu = 0
    # Z is independent of R, we effectively condition on Z
    # Unit vector used for finding the truncation interval
          # Fix scale
          # What is conditional variance?
    scale <- 1
    df <- sum(diag(P))
    # The truncated chi
    TC <- sqrt(sum(R^2))
    
    # For each step...
    L <- check_inequalities(fit, x, y, index, k, R)
    
    # Compute intersection:
    E <- do.call(interval_intersection, unlist(L, recursive = F))
    
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
    pvals <- c(pvals, value)    
  }
  names(pvals) <- fit$variable
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
      integral <- numerically_integrate(a, b, df)
  }
  return(integral)
}

numerically_integrate <- function(a, b, df, nsamp = 10000) {
  grid <- seq(from=a, to=b, length.out=nsamp)
  integrand <- dchisq(grid, df)
  return((b-a)*mean(integrand))
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
#' @param fit fitted model, the result of \link{groupfs}
#' @param scale Known/estimated scaling parameter
#' @param k Multiplier of model size, use \code{k = 2} for AIC, \code{k = log(n)} for BIC, or \code{k = 2log(p)} for RIC.
#' @return Matrix with two columns, edf and aic of fit, with rows indexed by step
extractAIC.groupfs <- function(fit, scale, k = 2, ...) {

  if (class(fit) != "groupfs") stop("Incorrect object type of fit")

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
#' @param fit fitted model, the result of \link{groupfs}
#' @return Matrix with two columns, edf and aic of fit, with rows indexed by step
model_select <- function(fit, alpha = .1, ...) {

  if (class(fit) != "groupfs") stop("Incorrect object type of fit")

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
#' @param fit fitted model, the result of \link{groupfs}
plot.groupfs <- function(fit, ...) {
  xrange <- 1:attr(fit, "steps")
  plot.default(x = xrange,  y = fit$p.value, xlab = "Variable", ylab = "P-value", cex = .5, xaxt = "n", pch = 19, ylim = c(0,1))
  points.default(x = xrange, y = fit$log$chisq, cex = .5)
  varnames <- attr(fit, "varnames")[fit$variable]
  if (!is.null(varnames)) {
    axis(1, at = xrange, labels=FALSE)
    text(x = xrange, y=-.08, labels = varnames, srt = 45, pos = 1, xpd = TRUE, cex = .4)
  }
}


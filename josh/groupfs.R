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
  maxprojs <- list()

  # Store other information from each step
  path.info <- data.frame(imax=integer(), k=integer(), L=numeric(), RSS=numeric(), RSSdrop=numeric(), chisq=numeric(), projections=numeric())

  modelrank <- 1
  
  # Begin main loop
  for (step in 1:steps) {
      
    added <- add1.groupfs(x.update, y.update, index, labels, inactive, k, Sigma)
    
    # Group to be added
    imax <- added$imax
    inactive <- setdiff(inactive, imax)
    active <- union(active, imax)    
    inactive.inds <- which(!index %in% active)

    # Rank of group
    added$k <- ncol(added$maxproj)
    modelrank <- modelrank + added$k

    # Stop without adding if model has become saturated
    if (modelrank >= n) break

    # Regress added group out of y and inactive x
    P.imax <- added$maxproj %*% t(added$maxproj)
    y.update <- y.update - P.imax %*% y.update
    x.update[, inactive.inds] <- x.update[, inactive.inds] - P.imax %*% x.update[, inactive.inds]

    projections[[step]] <- added$projections
    maxprojs[[step]] <- added$maxproj

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
    added$maxproj <- NULL
    path.info <- rbind(path.info, data.frame(added))
    
    if (verbose) print(added)
  }

  # Create output object
  value <- list(variable=path.info$imax, L=path.info$L, projections = projections, maxprojs = maxprojs, log = path.info)
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

  # Compute singular vectors of projections
  # X_i %*% X_i^\dagger
  projections <- lapply(keys, function(i) {
    inds <- which(index == i)
    xi <- x[,inds]
    return(svd(xi)$u)
  })

  names(projections) <- keys

  # Compute sums of squares to determine which group is added
  # penalized by rank of group if k > 0
  terms <- lapply(keys, function(i) {
    Py <- t(projections[[i]]) %*% y
    sum(Py^2) - ifelse(k > 0, k * ncol(projections[[i]]), 0)
  })

  # Maximizer = group to be added
  terms.maxind <- which.max(terms)
  imax <- inactive[terms.maxind]
  keyind <- which(keys == imax)
  maxproj <- projections[[keyind]]
  projections[[keyind]] <- NULL
  
  L <- terms[[terms.maxind]]

  return(list(imax=imax, L=L, projections = projections, maxproj = maxproj))
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

flatten <- function(L) {
    if (is.list(L[[1]])) return(unlist(L, recursive=FALSE))
    return(L)
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


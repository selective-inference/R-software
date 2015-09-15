# Forward stepwise implementation with groups and truncated \chi p-values
# -----------------------------------------------------------

#' Forward stepwise
#'
#' @param x Matrix of predictors (n by p)
#' @param y Vector of outcomes (length n)
#' @param index Group membership indicator of length p
#' @param maxsteps Maximum number of steps for forward stepwise
#' @param normalize Should the design matrix be normalized?
#' @param sp Size penalty used when comparing groups of different sizes.
#' @return Index \code{imax} of added group, residualized data, truncated chi p-value
groupfs <- function(x, y, index, maxsteps, intercept = TRUE, normalize = TRUE, sp = 0, ...) UseMethod("groupfs")

groupfs.default <- function(x, y, index, maxsteps, intercept = TRUE, normalize = TRUE, sp = 0, verbose=FALSE, ...) {

  p <- ncol(x)
  n <- nrow(x)

  # Assume no groups if index not specified
  if (missing(index)) stop("Missing argument: index.")

  # Group labels
  labels <- unique(index)
  G <- length(labels)
  inactive <- labels
  active <- c()
  if (missing(maxsteps) || maxsteps >= min(n, G)) maxsteps <- min(n-1, G)
  
  # Initialize copies of data for loop
  by <- mean(y)
  y.update <- y
  if (intercept) y.update <- y - by
  y.last <- y.update  
  x.update <- x

  # Center and scale design matrix
  xscaled <- scale_groups(x.update, index, scale = normalize)
  xm <- xscaled$xm
  xs <- xscaled$xs
  x.update <- xscaled$x

  x.begin <- x.update
  y.begin <- y.update
  # Store all projections computed along the path
  projections <- list()
  maxprojs <- list()

  # Store other information from each step
  path.info <- data.frame(imax=integer(), k=integer(), L=numeric(), RSS=numeric(), RSSdrop=numeric(), chisq=numeric(), projections=numeric())

  modelrank <- 1
  
  # Begin main loop
  for (step in 1:maxsteps) {
      
    added <- add1.groupfs(x.update, y.update, index, labels, inactive, sp)
    
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
    added$RSS <- sum(y.update^2)
    scale.chisq <- 1

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
  value <- list(action=path.info$imax, L=path.info$L, projections = projections, maxprojs = maxprojs, log = path.info, index = index, y = y.begin, x = x.begin, bx = xm, sx = xs, intercept = intercept)
  class(value) <- "groupfs"
  attr(value, "labels") <- labels
  attr(value, "index") <- index
  attr(value, "maxsteps") <- maxsteps
  attr(value, "sp") <- sp
  if (is.null(attr(x, "varnames"))) {
    attr(value, "varnames") <- colnames(x)
  } else {
    attr(value, "varnames") <- attr(x, "varnames")
  }

  invisible(value)
}

# -----------------------------------------------------------

#' Add one group to the model
#'
#' @param x Design matrix
#' @param y Response vector
#' @param index Group membership indicator of length p
#' @param inactive Indices of inactive groups
#' @return Index \code{imax} of added group, residualized data, truncated chi p-value
add1.groupfs <- function(x, y, index, labels, inactive, sp, ...) {

  # Use characters to avoid issues where
  # list() populates NULL lists in the positions
  # of the active variables
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
  # penalized by rank of group if sp > 0
  terms <- lapply(keys, function(i) {
    Py <- t(projections[[i]]) %*% y
    val <- sum(Py^2)
    if (sp > 0) val <- val - sp * ncol(projections[[i]])
    return(val)
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
scale_groups <- function(x, index, center = TRUE, scale = TRUE, ...) {
  keys <- unique(index)
  xm <- rep(0, ncol(x))
  xs <- rep(1, ncol(x))
  
  for (j in keys) {
    inds <- which(index == j)
    if (center) {
        xmj <- mean(x[, inds])
        xm[inds] <- xmj
        x[, inds] <- x[, inds] - xmj
    }
    if (scale) {
      normsq <- sum(x[, inds]^2)
      xsj <- sqrt(normsq)
      xs[inds] <- xsj
      if (xsj > 0) x[, inds] <- x[, inds] / xsj
    }

  }
  return(list(x=x, xm=xm, xs=xs))
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
    models$last <- fit$action[1:ind]
  }

  # ForwardStop
  renyi.transform <- - cumsum(log(1-fit$p.value)) / 1:length(fit$p.value)
  ind <- max_which(renyi.transform <= alpha)
  if (is.null(ind)) {
    models$forward <- 0
  } else {
    models$forward <- fit$action[1:ind]
  }

  # RIC
  ind <- which.min(extractAIC(fit, k = 2*log(p))[, 2])
  models$RIC <- fit$action[1:ind]

  # BIC
  ind <- which.min(extractAIC(fit, k = log(n))[, 2])
  models$BIC <- fit$action[1:ind]

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
  varnames <- attr(fit, "varnames")[fit$action]
  if (!is.null(varnames)) {
    axis(1, at = xrange, labels=FALSE)
    text(x = xrange, y=-.08, labels = varnames, srt = 45, pos = 1, xpd = TRUE, cex = .4)
  }
}


#' Select a model with forward stepwise.
#'
#' This function implements forward selection of linear models almost identically to \code{\link{stepAIC}} with \code{direction = "forward"}. The reason this is a separate function from \code{\link{fs}} is that groups of variables (e.g. dummies encoding levels of a categorical variable) must be handled differently in the selective inference framework.
#'
#' @param x Matrix of predictors (n by p).
#' @param y Vector of outcomes (length n).
#' @param index Group membership indicator of length p.
#' @param maxsteps Maximum number of steps for forward stepwise.
#' @param sigma Estimate of error standard deviation for use in AIC criterion. This determines the relative scale between RSS and the degrees of freedom penalty. Default is NULL corresponding to unknown sigma. See \code{\link{extractAIC}} for details.
#' @param k Multiplier of model size penalty, the default is \code{k = 2} for AIC. Use \code{k = log(n)} for BIC, or \code{k = log(p)} for RIC.
#' @param intercept Should an intercept be included in the model? Default is TRUE.
#' @param normalize Should the design matrix be normalized? Default is TRUE.
#' @param verbose Print out progress along the way? Default is FALSE.
#' @return An object of class "groupfs" containing information about the sequence of models in the forward stepwise algorithm. Call the function \code{\link{groupfsInf}} on this object to compute selective p-values.
#' @examples
#' x = matrix(rnorm(20*40), nrow=20)
#' index = sort(rep(1:20, 2))
#' y = rnorm(20) + 2 * (x[,1] - x[,2]) - (x[,3] - x[,4])
#' fit = groupfs(x, y, index, maxsteps = 5)
#' pvals = groupfsInf(fit)
#' @seealso \code{\link{groupfsInf}}
groupfs <- function(x, y, index, maxsteps, sigma = NULL, k = 2, intercept = TRUE, normalize = TRUE, verbose = FALSE) {

  if (missing(index)) stop("Missing argument: index.")
  p <- ncol(x)
  n <- nrow(x)

  # Group labels
  labels <- unique(index)
  G <- length(labels)
  inactive <- labels
  active <- c()

  if (missing(maxsteps) || maxsteps >= min(n, G)) maxsteps <- min(n-1, G)
  checkargs.xy(x=x, y=y)
  checkargs.groupfs(x, index, maxsteps)
  if (maxsteps > G) stop("maxsteps is larger than number of groups")
  gsizes <- sort(rle(sort(index))$lengths, decreasing = TRUE)
  if (sum(gsizes[1:maxsteps]) >= nrow(x)) {
      maxsteps <- max(which(cumsum(gsizes) < nrow(x)))
      warning(paste("If the largest groups are included the model will be saturated/overdetermined. To prevent this maxsteps has been changed to", maxsteps))
  }


  # Initialize copies of data for loop
  by <- mean(y)
  y.update <- y
  if (intercept) y.update <- y - by
  y.last <- y.update
  x.update <- x

  # Center and scale design matrix
  xscaled <- scaleGroups(x.update, index, scale = normalize)
  xm <- xscaled$xm
  xs <- xscaled$xs
  x.update <- xscaled$x

  x.begin <- x.update
  y.begin <- y.update
  # Store all projections computed along the path
  terms = projections = maxprojs = aicpens = maxpens = cumprojs = vector("list", maxsteps)

  # Store other information from each step
  path.info <- data.frame(imax=integer(maxsteps), L=numeric(maxsteps), df=integer(maxsteps), RSS=numeric(maxsteps), RSSdrop=numeric(maxsteps), chisq=numeric(maxsteps))

  modelrank <- 1

  # Begin main loop
  for (step in 1:maxsteps) {

    added <- add1.groupfs(x.update, y.update, index, labels, inactive, k, sigma)

    # Group to be added
    imax <- added$imax
    inactive <- setdiff(inactive, imax)
    active <- union(active, imax)
    inactive.inds <- which(!index %in% active)

    # Rank of group
    modelrank <- modelrank + added$df

    # Stop without adding if model has become saturated
    if (modelrank >= n) {
        stop("Saturated model. Abandon ship!")
    }

    # Regress added group out of y and inactive x
    P.imax <- added$maxproj %*% t(added$maxproj)
    if (is.null(sigma)) {
        P.imax <- P.imax / exp(k*added$df/n)
    } else {
        P.imax <- P.imax * sigma^2
    }
    P.imax <- diag(rep(1, n)) - P.imax
    y.update <- P.imax %*% y.update
    x.update[, inactive.inds] <- P.imax %*% x.update[, inactive.inds]

    projections[[step]] <- added$projections
    maxprojs[[step]] <- added$maxproj
    aicpens[[step]] <- added$aicpens
    maxpens[[step]] <- added$maxpen
    if (step == 1) cumprojs[[step]] <- P.imax
    if (step > 1) cumprojs[[step]] <- P.imax %*% cumprojs[[step-1]]
    terms[[step]] <- added$terms

    # Compute RSS for unadjusted chisq p-values
    added$RSS <- sum(y.update^2)
    scale.chisq <- 1

    added$RSSdrop <- sum((y.last - y.update)^2)
    added$chisq <- pchisq(added$RSSdrop/scale.chisq, lower.tail=FALSE, df = added$df)
    y.last <- y.update

    # Projections are stored separately
    step.info <- data.frame(added[-c(4:(length(added)-3))])
    path.info[step, ] <- step.info

    if (verbose) print(step.info)
  }

  # Create output object
  value <- list(action=path.info$imax, L=path.info$L, projections = projections, maxprojs = maxprojs, aicpens = aicpens, maxpens = maxpens, cumprojs = cumprojs, log = path.info, index = index, y = y.begin, x = x.begin, bx = xm, sx = xs, sigma = sigma, intercept = intercept, terms = terms)
  class(value) <- "groupfs"
  attr(value, "labels") <- labels
  attr(value, "index") <- index
  attr(value, "maxsteps") <- maxsteps
  attr(value, "sigma") <- sigma
  attr(value, "k") <- k
  if (is.null(attr(x, "varnames"))) {
    attr(value, "varnames") <- colnames(x)
  } else {
    attr(value, "varnames") <- attr(x, "varnames")
  }

  invisible(value)
}

#' Add one group to the model in \code{groupfs}.
#'
#' For internal use by \code{\link{groupfs}}.
#'
#' @param x Design matrix.
#' @param y Response vector.
#' @param index Group membership indicator of length p.
#' @param labels The unique elements of \code{index}.
#' @param inactive Labels of inactive groups.
#' @param k Multiplier of model size penalty, use \code{k = 2} for AIC, \code{k = log(n)} for BIC, or \code{k = log(p)} for RIC.
#' @param sigma Estimate of error standard deviation for use in AIC criterion. This determines the relative scale between RSS and the degrees of freedom penalty. See \code{\link{extractAIC}} for details.
#' @return Index \code{imax} of added group, value \code{L} of maximized negative AIC, lists of projection matrices defining quadratic model selection event.
add1.groupfs <- function(x, y, index, labels, inactive, k, sigma = NULL) {

  # Use characters to avoid issues where
  # list() populates NULL lists in the positions
  # of the active variables
  ### Question for later: does this slow down lapply?
  keys = as.character(inactive)
  n2y <- sum(y^2)
  n <- ncol(x)

  # Compute sums of squares to determine which group is added
  # penalized by rank of group if k > 0
  projections = aicpens = terms = vector("list", length(keys))
  names(projections) = names(terms) = names(aicpens) = keys
  for (key in keys) {
      inds <- which(index == key)
      xi <- x[,inds]
      ui <- svdu_thresh(xi)
      dfi <- ncol(ui)
      projections[[key]] <- ui
      dfi <- ncol(ui)
      uy <- t(ui) %*% y
      if (is.null(sigma)) {
          aicpens[[key]] <- exp(k*dfi/n)
          terms[[key]] <- (sum(uy^2)  - sum(y^2)) * aicpens[[key]]
      } else {
          aicpens[[key]] <- sigma^2 * k * dfi/n
          terms[[key]] <- (sum(uy^2) - sum(y^2)) - aicpens[[key]]
      }
  }

  # Maximizer = group to be added
  terms.maxind <- which.max(terms)
  imax <- inactive[terms.maxind]
  keyind <- which(keys == imax)
  maxproj <- projections[[keyind]]
  maxpen <- aicpens[[keyind]]
  projections[[keyind]] <- NULL
  aicpens[[keyind]] <- NULL

  L <- terms[[terms.maxind]]

  return(list(imax=imax, L=L, df = ncol(maxproj), projections = projections, maxproj = maxproj, aicpens = aicpens, maxpen = maxpen, terms = terms))
}

# -----------------------------------------------------------

#' Compute selective p-values for a model fitted by \code{groupfs}.
#'
#' Computes p-values for each group of variables in a model fitted by \code{\link{groupfs}}. These p-values adjust for selection by truncating the usual \code{chi^2} statistics to the regions implied by the model selection event. Details are provided in a forthcoming work.
#'
#' @param obj Object returned by \code{\link{groupfs}} function
#' @param sigma Estimate of error standard deviation. If NULL (default), this is estimated using the mean squared residual of the full least squares fit when n >= 2p, and the mean squared residual of the selected model when n < 2p. In the latter case, the user should use \code{\link{estimateSigma}} function for a more accurate estimate.
#' @param projs Additional projections to define model selection event. For use with cross-validation. Default is NULL and it is not recommended to change this.
#' @param verbose Print out progress along the way? Default is FALSE.
#' @return An object of class "groupfsInf" containing selective p-values for the fitted model \code{obj}. The default printing behavior should supply adequate information.
#'
#' \describe{
#'   \item{vars}{Labels of the active groups in the order they were included.}
#'   \item{pv}{Selective p-values computed from appropriate truncated distributions.}
#'   \item{sigma}{Estimate of error variance used in computing p-values.}
#'   \item{TC}{Observed value of truncated chi.}
#'   \item{df}{Rank of group of variables when it was added to the model.}
#'   \item{support}{List of intervals defining the truncation region of the truncated chi.}
#' }
groupfsInf <- function(obj, sigma = NULL, projs = NULL, verbose = FALSE) {

  n <- nrow(obj$x)
  p <- ncol(obj$x)
  active <- obj$action
  maxsteps <- attr(obj, "maxsteps")
  k <- attr(obj, "k")
  index <- obj$index
  x <- obj$x
  y <- obj$y
  Eindex <- which(index %in% active)
  Ep <- length(Eindex)

  nanconv <- FALSE
  pvals <- numeric(maxsteps)
  dfs <- numeric(maxsteps)
  TCs <- numeric(maxsteps)
  supports <- list()

  if (!is.null(sigma)) {
      if (!is.null(obj$sigma)) {
          cat(paste("Using specified value", sigma, "for sigma in place of the value", obj$sigma, "used by groupfs()\n"))
      }
  } else {
      if (is.null(obj$sigma)) {
          if (n >= 2*p) {
              sigma <- sqrt(sum(lsfit(obj$x, obj$y, intercept = obj$intercept)$res^2)/(n-p-obj$intercept))
          } else {
              sigma = sqrt(obj$log$RSS[length(obj$log$RSS)]/(n-Ep-obj$intercept))
              warning(paste(sprintf("p > n/2, and sigmahat = %0.3f used as an estimate of sigma;",sigma), "you may want to use the estimateSigma function"))
          }
      } else {
          sigma <- obj$sigma
      }
  }

  # Compute p-value for each active group
  for (j in 1:maxsteps) {
    i <- active[j]
    if (verbose) cat(paste0("Step ", j, "/", maxsteps, ": computing P-value for group ", i, "\n"))    
    # Form projection onto active set minus i
    # and project x_i orthogonally
    x_i <- x[,which(index == i), drop = FALSE]
    if (length(active) > 1) {
        minus_i <- setdiff(active, i)
        x_minus_i <- svdu_thresh(x[,which(index %in% minus_i), drop = FALSE])
        x_i <- x_i - x_minus_i %*% t(x_minus_i) %*% x_i
    }

    # Project y onto what remains of x_i
    Ugtilde <- svdu_thresh(x_i)
    R <- t(Ugtilde) %*% y
    TC <- sqrt(sum(R^2))
    eta <- Ugtilde %*% R / TC
    df <- ncol(Ugtilde)
    dfs[j] <- df
    TCs[j] <- TC

    # For each step...
    L <- interval_groupfs(obj, TC, R, eta, Ugtilde)

    # Any additional projections, e.g. from cross-validation?
    if (!is.null(projs)) L <- c(L, projs)

    # Compute intersection:
    Lunion <- do.call(interval_union, L)
    Lunion <- interval_union(Lunion, Intervals(c(-Inf,0)))
    E <- interval_complement(Lunion, check_valid = FALSE)
    supports[[j]] <- E

    # E is now potentially a union of intervals
    if (length(E) == 0) stop("Trivial intersection")

    # Sum truncated cdf over each part of E
    denom <- do.call(sum, lapply(1:nrow(E), function(v) {
      tchi_interval(E[v,1], E[v,2], sigma, df)
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
          return(tchi_interval(TC, upper, sigma, df))
        } else {
          # Observed value is not in this interval
          return(tchi_interval(lower, upper, sigma, df))
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
      nanconv <- TRUE
    }
    # Force p-value to lie in the [0,1] interval
    # in case of numerical issues
    value <- max(0, min(1, value))
    pvals[j] <- value
  }
  if (nanconv) warning("P-value NaNs of the form 0/0 converted to 0. This typically occurs for numerical reasons in the presence of a large signal-to-noise ratio.")
  names(pvals) <- obj$action
  out <- list(vars = active, pv=pvals, sigma=sigma, TC=TCs, df = dfs, support=supports)
  class(out) <- "groupfsInf"
  if (!is.null(attr(obj, "varnames"))) {
      attr(out, "varnames") <- attr(obj, "varnames")
  }
  invisible(out)
}

# -----------------------------------------------------------

tchi_interval <- function(lower, upper, sigma, df) {
  a <- (lower/sigma)^2
  b <- (upper/sigma)^2
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


#' Center and scale design matrix by groups
#'
#' @param x Design matrix.
#' @param index Group membership indicator of length p.
#' @param center Center groups, default is TRUE.
#' @param scale Scale groups by Frobenius norm, default is TRUE.
#' @return Scaled design matrix
scaleGroups <- function(x, index, center = TRUE, scale = TRUE) {
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
    normsq <- sum(x[, inds]^2)
    xsj <- sqrt(normsq)
    xs[inds] <- xsj
    if (xsj > 0) {
        if (scale) x[, inds] <- x[, inds] / xsj
    } else {
        stop(paste("Design matrix contains identically zero group of variables:", j))
    }
  }
  return(list(x=x, xm=xm, xs=xs))
}

#' Expand a data frame with factors to form a design matrix with the full binary encoding of each factor.
#' 
#' When using \code{\link{groupfs}} with factor variables call this function first to create a design matrix. 
#'
#' @param df Data frame containing some columns which are \code{factors}.
#' @return List containing
#' \describe{
#'   \item{x}{Design matrix, the first columns contain any numeric variables from the original date frame.}
#'   \item{index}{Group membership indicator for expanded matrix.}
#' }
#' @examples
#' \dontrun{
#' fd = factorDesign(warpbreaks)
#' y = rnorm(nrow(fd$x))
#' fit = groupfs(fd$x, y, fd$index, maxsteps=2, intercept=F)
#' pvals = groupfsInf(fit)
#' }
factorDesign <- function(df) {
    factor.inds <- sapply(df[1,], is.factor)
    factor.labels <- which(factor.inds)
    nfacs <- sum(factor.inds)
    nlevs <- sapply(df[1,factor.inds], function(fac) nlevels(fac))
    totnlevs <- sum(nlevs)
    num.num = indcounter = ncol(df) - nfacs
    x <- matrix(nrow=nrow(df), ncol = totnlevs + num.num)
    colnames(x) <- 1:ncol(x)
    index <- integer(ncol(x))
    varnames <- character(ncol(df))
    if (num.num > 0) {
        x[,1:num.num] <- df[, !factor.inds]
        varnames[1:num.num] = colnames(x)[1:num.num] <- colnames(df)[1:num.num]
        index[1:num.num] <- 1:num.num
        indcounter <- indcounter + num.num - 1
    }
    for (j in 1:nfacs) {
        submat <- model.matrix(~ df[, factor.labels[j]] - 1)
        indcounter <- indcounter+1
        submatinds <- indcounter:(indcounter+nlevs[j]-1)
        indcounter <- indcounter + nlevs[j] - 1
        colnames(x)[submatinds] <- paste0(colnames(df)[num.num + j], ":", 1:nlevs[j])
        varnames[num.num + j] <- colnames(df)[num.num + j]
        x[,submatinds] <- submat
        index[submatinds] <- num.num + j
    }
    attr(x, "varnames") <- varnames
    return(list(x = x, index = index))
}

svdu_thresh <- function(x) {
    svdx <- svd(x)
    inds <- svdx$d > svdx$d[1] * sqrt(.Machine$double.eps)
    return(svdx$u[, inds, drop = FALSE])
}

flatten <- function(L) {
    if (is.list(L[[1]])) return(unlist(L, recursive=FALSE))
    return(L)
}

print.groupfs <- function(x, ...) {
    cat("\nSequence of added groups:\n")
    nsteps = length(x$action)
    action <- x$action
    vnames <- attr(x, "varnames")
    if (length(vnames) > 0) action <- vnames[action]
    tab = data.frame(Group = action, Rank = x$log$df, RSS = round(x$log$RSS, 3))
    rownames(tab) = 1:nsteps
    print(tab)
    cat("\nUse groupfsInf() to compute P-values\n")
    invisible()
}

print.groupfsInf <- function(x, ...) {
    cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n", x$sigma))
    action <- x$vars
    vnames <- attr(x, "varnames")
    if (length(vnames) > 0) action <- vnames[action]
    tab = data.frame(Group = action, Pvalue = round(x$pv, 3), Tchi = round(x$TC, 3),
        df = x$df, Size = round(unlist(lapply(lapply(x$support, size), sum)), 3),
        Ints = unlist(lapply(x$support, nrow)), Min =round(unlist(lapply(x$support, min)), 3),
        Max = round(unlist(lapply(x$support, max)), 3))
    rownames(tab) = 1:length(x$vars)
    print(tab)
    cat("\nInts is the number of intervals in the truncated chi selection region and Size is the sum of their lengths. Min and Max are the lowest and highest endpoints of the truncation region. No confidence intervals are reported by groupfsInf.\n")
    invisible()
}

checkargs.groupfs <- function(x, index, maxsteps) {
    if (length(index) != ncol(x)) stop("Length of index does not match number of columns of x")
    if ((round(maxsteps) != maxsteps) || (maxsteps <= 0)) stop("maxsteps must be an integer > 0")
}

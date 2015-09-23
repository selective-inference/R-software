#' Select a model with forward stepwise.
#'
#' \code{groupfs} implements forward selection of linear models almost identically to \code{\link{stepAIC}} with \code{direction = "forward"}. The reason this is a separate function from \code{\link{fs}} is that groups of variables (e.g. dummies encoding levels of a categorical variable) must be handled differently in the selective inference framework.
#'
#' @param x Matrix of predictors (n by p).
#' @param y Vector of outcomes (length n).
#' @param index Group membership indicator of length p.
#' @param maxsteps Maximum number of steps for forward stepwise.
#' @param sigma Estimate of error standard deviation for use in AIC criterion. This determines the relative scale between RSS and the degrees of freedom penalty. Default is NULL corresponding to unknown sigma. See \code{\link{extractAIC}} for details.
#' @param k Multiplier of model size penalty, use \code{k = 2} for AIC, \code{k = log(n)} for BIC, or \code{k = log(p)} for RIC.
#' @param intercept Should an intercept be included in the model? Default is TRUE.
#' @param normalize Should the design matrix be normalized? Default is TRUE.
#' @param verbose Print out progress along the way? Default is FALSE.
#' @return An object of class "groupfs" containing information about the sequence of models in the forward stepwise algorithm. Call the function \code{\link{groupfsInf}} on this object to compute selective p-values.
#' @examples
#' x = matrix(rnorm(10*40, nrow=10))
#' index = sort(rep(1:20, 2))
#' y = rnorm(10) + x[,1] - x[,2]
#' fit = groupfs(x, y, index)
#' pvals = groupfsInf(fit)
#' @seealso \code{\link{groupfsInf}}
groupfs <- function(x, y, index, maxsteps, sigma = NULL, k = 2, intercept = TRUE, normalize = TRUE, verbose = FALSE) UseMethod("groupfs")

groupfs.default <- function(x, y, index, maxsteps, sigma = NULL, k = 2, intercept = TRUE, normalize = TRUE, verbose=FALSE) {

  p <- ncol(x)
  n <- nrow(x)
  checkargs.xy(x=x, y=y)
  checkargs.groupfs(x, index, maxsteps)

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
  projections = maxprojs = aicpens = maxpens = cumprojs = vector("list", maxsteps)

  # Store other information from each step
  path.info <- data.frame(imax=integer(maxsteps), L=numeric(maxsteps), df=integer(maxsteps), RSS=numeric(maxsteps), RSSdrop=numeric(maxsteps), chisq=numeric(maxsteps))

  modelrank <- 1
  cumprojs[[1]] <- diag(rep(1, n))

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
    if (modelrank >= n) break

    # Regress added group out of y and inactive x
    P.imax <- added$maxproj %*% t(added$maxproj)
    if (is.null(sigma)) {
        P.imax <- P.imax / exp(k*added$df/n)
    } else {
        P.imax <- P.imax * sigma^2
    }
    y.update <- y.update - P.imax %*% y.update
    x.update[, inactive.inds] <- x.update[, inactive.inds] - P.imax %*% x.update[, inactive.inds]

    projections[[step]] <- added$projections
    maxprojs[[step]] <- added$maxproj
    aicpens[[step]] <- added$aicpens
    maxpens[[step]] <- added$maxpen
    if (step == 1) cumprojs[[step]] <- P.imax
    if (step > 1) cumprojs[[step]] <- cumprojs[[step-1]] %*% P.imax

    # Compute RSS for unadjusted chisq p-values
    added$RSS <- sum(y.update^2)
    scale.chisq <- 1

    added$RSSdrop <- sum((y.last - y.update)^2)
    added$chisq <- pchisq(added$RSSdrop/scale.chisq, lower.tail=FALSE, df = added$df)
    y.last <- y.update

    # Projections are stored separately
    step.info <- data.frame(added[-c(4:7)])
    path.info[step, ] <- step.info

    if (verbose) print(step.info)
  }

  # Create output object
  value <- list(action=path.info$imax, L=path.info$L, projections = projections, maxprojs = maxprojs, aicpens = aicpens, maxpens = maxpens, cumprojs = cumprojs, log = path.info, index = index, y = y.begin, x = x.begin, bx = xm, sx = xs, intercept = intercept)
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

# -----------------------------------------------------------

#' Add one group to the model. For internal use by \code{\link{groupfs}}.
#'
#' @param x Design matrix.
#' @param y Response vector.
#' @param index Group membership indicator of length p.
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

  # Compute singular vectors of projections
  # X_i %*% X_i^\dagger
  projections <- lapply(keys, function(i) {
    inds <- which(index == i)
    xi <- x[,inds]
    svdi <- svd(xi)
    #inds <- svdi$d > svdi$d[1] * sqrt(.Machine$double.eps)
    inds <- 1:ncol(svdi$u)
    ui <- svdi$u[, inds, drop = FALSE]
    dfi <- ncol(ui)
    if (is.null(sigma)) {
        ui <- ui * exp(k*dfi/(2*n))
    } else {
        ui <- ui / sigma
    }
    return(ui)
  })

  names(projections) <- keys

  # Compute sums of squares to determine which group is added
  # penalized by rank of group if k > 0
  aicpens = terms = vector("list", length(keys))
  names(terms) = names(aicpens) = keys
  for (key in keys) {
      Ui <- projections[[key]]
      dfi <- ncol(Ui)
      Uy <- t(Ui) %*% y
      if (is.null(sigma)) {
          aicpens[[key]] <- n2y * exp(k*dfi/n)
          terms[[key]] <- sum(Uy^2)  - aicpens[[key]]
      } else {
          aicpens[[key]] <- k * dfi/n #- n2y / sigma^2
          terms[[key]] <- sum(Uy^2) - aicpens[[key]]
      }
  }

  # Maximizer = group to be added
  terms.maxind <- which.max(terms)
  imax <- inactive[terms.maxind]
  print(imax)
  keyind <- which(keys == imax)
  maxproj <- projections[[keyind]]
  maxpen <- aicpens[[keyind]]
  projections[[keyind]] <- NULL
  aicpens[[keyind]] <- NULL

  L <- terms[[terms.maxind]]

  return(list(imax=imax, L=L, df = ncol(maxproj), projections = projections, maxproj = maxproj, aicpens = aicpens, maxpen = maxpen))
}

# -----------------------------------------------------------

#' Form univariate selection interval for a given contrast
#'
#' @param obj Object returned by \code{\link{groupfs}} function
#' @param sigma Estimate of error standard deviation. If NULL (default), this is estimated using the mean squared residual of the full least squares fit when n >= 2p, and the mean squared residual of the selected model when n < 2p. In the latter case, the user should use \code{\link{estimateSigma}} function for a more accurate estimate.
#' @param projs Additional projections to define model selection event. For use with cross-validation. Default is NULL and it is not recommended to change this.
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
groupfsInf <- function(obj, sigma = NULL, projs = NULL) {

  cat(paste("\nUsing sigma value:", attr(obj, "sigma"), "\n"))

  x <- obj$x
  n <- nrow(x)
  p <- ncol(x)
  y <- obj$y
  active <- obj$action
  maxsteps <- attr(obj, "maxsteps")
  k <- attr(obj, "k")
  index <- obj$index
  Eindex <- which(index %in% active)
  Ep <- length(Eindex)

  nanconv <- FALSE
  pvals <- numeric(maxsteps)
  dfs <- numeric(maxsteps)
  TCs <- numeric(maxsteps)
  supports <- list()

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
    # Form projection onto active set minus i
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
    dfs[j] <- df
    # The truncated chi
    TC <- sqrt(sum(R^2))
    TCs[j] <- TC

    # For each step...
    L <- interval_groupfs(obj$action, obj$projections, obj$maxprojs, obj$cumprojs, x, y, index, k, obj$sigma, TC, R, Ugtilde)

    # Any additional projections, e.g. from cross-validation?
    if (!is.null(projs)) L <- c(L, projs)

    # Compute intersection:
    E <- interval_complement(do.call(interval_union, L), check_valid = FALSE)
    supports[[j]] <- E

    # E is now potentially a union of intervals
    if (length(E) == 0) stop("Trivial intersection")

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
  invisible(out)
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


#' Center and scale design matrix by groups
#'
#' @param x Design matrix.
#' @param index Group membership indicator of length p.
#' @param center Center groups, default is TRUE.
#' @param scale Scale groups by Frobenius norm, default is TRUE.
#' @return Scaled design matrix
scale_groups <- function(x, index, center = TRUE, scale = TRUE) {
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

print.groupfs <- function(x, ...) {
    cat("\nSequence of added groups:\n")
    nsteps = length(x$action)
    tab = cbind(1:nsteps, x$action, x$log$df)
    colnames(tab) = c("Step", "Group", "Rank")
    rownames(tab) = rep("", nrow(tab))
    print(tab)
    cat("\nUse groupfsInf() function to compute P-values\n")
    invisible()
}

print.groupfsInf <- function(x, ...) {
    cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
                              x$sigma))
    tab = cbind(x$vars, round(x$pv, 3), round(x$TC, 3), x$df,
        round(unlist(lapply(lapply(pvals$support, size), sum)), 3),
        unlist(lapply(pvals$support, nrow)), round(unlist(lapply(pvals$support, min)), 3),
        round(unlist(lapply(pvals$support, max)), 3))
    colnames(tab) = c("Var", "P-value", "Tchi", "df", "Int. size", "Components",
                "Int. inf", "Int. sup")
    rownames(tab) = rep("", nrow(tab))
    print(tab)
    cat("\nInt. inf and Int. sup are the lowest and highest endpoints of the truncation interval. No confidence intervals are reported by groupfsInf.\n")
    invisible()
}

checkargs.groupfs <- function(x, index, maxsteps) {
    if (missing(index)) stop("Missing argument: index.")
    if (length(index) != ncol(x)) stop("Length of index does not match number of columns of x")
    G <- length(unique(index))
    if (maxsteps > G) stop("maxsteps is larger than number of groups")
    gsizes <- sort(rle(sort(index))$lengths, decreasing = TRUE)
    if (sum(gsizes[1:maxsteps]) > nrow(x)) stop("maxsteps is too large. If the largest groups are included the model is overdetermined")
}


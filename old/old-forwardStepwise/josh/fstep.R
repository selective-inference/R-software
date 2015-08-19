# Forward stepwise implementation with groups, weights, and Kac-Rice p-values

# -----------------------------------------------------------

#' Forward stepwise, computing tchi p-values along the path
#'
#' @param x Design matrix
#' @param y Response vector
#' @param index Group membership indicator of length p
#' @param wt Pre-specified weight for each group
#' @param Sigma Error covariance matrix
#' @param steps Maximum number of steps for forward stepwise
#' @param normalize Should the design matrix be normalized?
#' @param mc.cores Number of cores for parallel computing
#' @return Index \code{imax} of added group, residualized data, truncated chi p-value
fstep <- function(x, y, index, wt, Sigma, steps, normalize = TRUE, mc.cores = 1, ...) UseMethod("fstep")

fstep.default <- function(x, y, index, wt, Sigma, steps, normalize = TRUE, mc.cores = 1, verbose=FALSE, ...) {

  p <- ncol(x)
  n <- nrow(x)

  if (missing(index)) index <- 1:p
  labels <- unique(index)
  if (missing(wt)) wt <- c(rep(1, length(labels)))
  if (missing(Sigma)) {
    Sigma <- 1
    warning("Error covariance unspecified, using diag(1)")
  }
  if (missing(steps)) steps <- min(n, p) - 1
  inactive <- labels

  # Form whitening matrix Sigma^(-1/2)
  if (is.matrix(Sigma)) {
    svdSigma <- svd(Sigma)
    whitener <- svdSigma$u %*% diag(1/sqrt(svdSigma$d)) %*% t(svdSigma$u)
  }
  active <- c()
  x.update <- x
  if (normalize) {
    x.update <- scale_groups(x.update, index)
  }

  y.update <- y - mean(y)
  y.last <- y.update

  residmat <- matrix(NA, nrow = steps, ncol = length(y))
  x.active <- matrix(NA, nrow = n, ncol = 1)

  output <- data.frame(imax=integer(), p.value=numeric(), k=integer(), conditional_var=numeric(), L=numeric(), lower_bound=numeric(), upper_bound=numeric(), RSS=numeric(), RSSdrop=numeric(), chisq=numeric())

  for (step in 1:steps) {
    added <- add1.fstep(x.update, y.update, index, inactive, wt, Sigma) #, ...)
    imax <- added$imax
    inactive <- setdiff(inactive, imax)
    active.inds <- which(index %in% active)
    inactive.inds <- which(!index %in% active)
    imax.inds <- which(index == imax)
    ## if (length(active) > 0) {
    ##   x.imax <- lm(x.update[, imax.inds] ~ x.active[, -1] - 1)$residuals
    ## } else {
      x.imax <- x.update[, imax.inds]
    ## }
    y.update <- resid(lm(y.update ~ x.imax - 1))
    x.update[, inactive.inds] <- resid(lm(x.update[, inactive.inds] ~ x.imax - 1))
    ## xginv <- ginv(x.imax)
    ## P.imax <- x.imax %*% xginv
    ## y.update <- y.update - P.imax %*% y.update
    ## x.update[, inactive.inds] <- x.update[, inactive.inds] - P.imax %*% x.update[, inactive.inds]

    active <- union(active, imax)
    x.active <- cbind(x.active, x.imax)
    residmat[step, ] <- y.update

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
    output <- rbind(output, data.frame(added))
    if (verbose) print(added)
  }

  value <- list(variable=output$imax, p.value=output$p.value, log=output)
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
#' @param wt Pre-specified weight for each group
#' @param Sigma Error covariance matrix
#' @param mc.cores Number of cores for parallel computing
#' @return Index \code{imax} of added group, residualized data, truncated chi p-value
add1.fstep <- function(x, y, index, inactive, wt, Sigma, mc.cores = 1, ...) {

  n <- nrow(x)
  p <- ncol(x)
  Xy <- t(x) %*% y

  if (missing(index)) index <- 1:p
  labels <- unique(index)
  if (missing(wt)) wt <- c(rep(1, length(labels)))
  if (missing(inactive)) inactive <- labels
  if (missing(Sigma)) {
    Sigma <- 1
    warning("Error covariance unspecified, using diag(1)")
  }

  wt.list <- list()
  wt.list[labels] <- wt

  terms <- lapply(inactive, function(i) {
    inds <- which(index == i)
    return(sum((Xy[inds])^2)/wt.list[[i]]^2)
  })

  terms.maxind <- which.max(terms)
  imax <- inactive[terms.maxind]
  maxinds <- which(index == imax)
  L <- sqrt(terms[[terms.maxind]])
  w_max <- wt.list[[imax]]
  Xy_max <- Xy[maxinds]
  Xmax <- x[, maxinds]

  active.inds <- which(!index %in% inactive)
  inactive.inds <- which(index %in% setdiff(inactive, imax))

  Xmax_ncol <- ncol(Xmax)
  if (is.null(Xmax_ncol)) Xmax_ncol <- 1

  if (is.matrix(Xmax)) {
    k <- rankMatrix(Xmax)[1]
  } else {
    k <- 1
  }

  eta <- rep(0, p)
  eta[maxinds] <- Xy_max / (sqrt(sum(Xy_max^2)) * w_max)

  if (Xmax_ncol > 1) {
    # Group of size > 1, non-trivial gradient term
    tangent.space <- Null(Xy_max)
    tangent.space = cbind(tangent.space, rep(0, Xmax_ncol))
    V <- matrix(0, ncol = Xmax_ncol, nrow = p)
    V[maxinds, ] <- tangent.space
    XV <- x %*% V[, -ncol(V)]
    XV <- cbind(XV, rep(0, nrow(XV)))
    XXy <- Xmax %*% Xy_max

    if (is.matrix(Sigma)) {
      SXV <- Sigma %*% XV
    } else {
      SXV <- Sigma*XV
    }

    VXSXV <- t(XV) %*% SXV
    #Categorical designs and large signals cause trouble here
    if (max(abs(VXSXV)) < .Machine$double.eps) {
      OP <- 1
    } else {
       P <- SXV %*% ginv(VXSXV) %*% t(XV)
       OP <- diag(rep(1, n)) - P
    }

    if ((is.matrix(Sigma)) && (is.matrix(OP))) {
      H <- OP %*% Sigma
    } else {
      H <- OP*Sigma
    }

    if (is.matrix(H)) {
      HXXy <- H %*% XXy
    } else {
      HXXy <- H*XXy
    }

    temp_den <- w_max * sqrt(sum(Xy_max^2))
    Xeta <- HXXy / temp_den
    conditional_var <- t(XXy) %*% HXXy
    # diag() coerces 1x1 matrix to numeric
    conditional_var <- diag(conditional_var) / temp_den^2

  } else {
    # Group of size 1, no gradient term
    Xeta <- Xmax / w_max * sign(Xy_max)
    if (is.matrix(Sigma)) {
      conditional_var <- t(Xeta) %*% Sigma %*% Xeta
    } else {
      conditional_var <- Sigma * sum(Xeta^2)
    }
  }

  if (conditional_var <= 0) {
    stop(paste("Conditional variance negative:", conditional_var))
    p.value <- NaN
    lower_bound <- NaN
    upper_bound <- NaN

  } else {
    Xeta <- Xeta / conditional_var
    b <- t(x) %*% Xeta
    a <- Xy - b * L

    inactive_new <- setdiff(inactive, imax)

    if (mc.cores > 1) {
      trig_solutions <- mclapply(inactive_new, function(i) {
        inds <- which(index == i)
        w <- wt.list[[i]]
        return(tchi_trig_solution(num=a[inds], den=b[inds], w=w))
      }, mc.cores = mc.cores)

    } else {
      trig_solutions <- lapply(inactive_new, function(i) {
        inds <- which(index == i)
        w <- wt.list[[i]]
        return(tchi_trig_solution(num=a[inds], den=b[inds], w=w))
      })
    }

    trig_solutions <- do.call(cbind, trig_solutions)
    V_lower <- trig_solutions[1, ]
    V_upper <- trig_solutions[2, ]

    if (length(V_lower) >= 1) {
      V_lower <- V_lower[!is.nan(V_lower)]
      lower_bound <- max(V_lower)
    } else {
      lower_bound <- 0
    }
    if (length(V_upper) >= 1) {
      V_upper <- V_upper[!is.nan(V_upper)]
      upper_bound <- min(V_upper)
    } else {
      upper_bound <- Inf
    }

    p.value <- tchi_cdf(L, lower_bound, upper_bound, conditional_var, k)
  }

  return(list(imax=imax, p.value=p.value, k=k, conditional_var=conditional_var, L=L, lower_bound=lower_bound, upper_bound=upper_bound))
}

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

  # Avoid infuriating integer(0) issues with which()
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


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
#' @return Index \code{imax} of added group, residualized data, truncated chi p-value
fstep <- function(x, y, index, wt, Sigma, steps, normalize = TRUE, ...) UseMethod("fstep")

fstep.default <- function(x, y, index, wt, Sigma, steps, normalize = TRUE, verbose=FALSE, ...) {

  p <- ncol(x)
  n <- nrow(x)

  if (missing(index)) index <- 1:p
  labels <- unique(index)
  if (missing(wt)) wt <- c(rep(1, length(labels)))
  if (missing(Sigma)) {
    Sigma <- 1
    if (verbose) warning("Error covariance unspecified, using diag(1)")
  }
  if (missing(steps)) steps <- min(n, p) - 1
  wt.list <- list()
  wt.list[labels] <- wt  
  inactive <- labels

  # Form whitening matrix Sigma^(-1/2)
  if (is.matrix(Sigma)) {
    svdSigma <- svd(Sigma)
    whitener <- svdSigma$u %*% diag(1/sqrt(svdSigma$d)) %*% t(svdSigma$u)
  }
  active <- c()
  xr <- x
  if (normalize) {
    xr <- scale_groups(xr, index)
  }

  r <- y - mean(y)
  r.last <- r

  output <- data.frame(imax=integer(), p.value=numeric(), k=integer(), conditional_var=numeric(), L=numeric(), lower_bound=numeric(), upper_bound=numeric(), RSS=numeric(), RSSdrop=numeric(), chisq=numeric())

  rs <- etas <- matrix(NA, nrow=n, ncol=steps)
  qs <- list()
  
  for (step in 1:steps) {

    rs[, step] <- r
    Xy <- t(xr) %*% r
    terms <- lapply(inactive, function(i) {
      inds <- which(index == i)
      return(sum((Xy[inds])^2)/wt.list[[i]]^2)
    })

    terms.maxind <- which.max(terms)
    imax <- inactive[terms.maxind]
    maxinds <- which(index == imax)
    #L <- sqrt(terms[[terms.maxind]])
    w_max <- wt.list[[imax]]

    inactive <- setdiff(inactive, imax)
    active.inds <- which(index %in% active)
    inactive.inds <- which(!index %in% active)
    imax.inds <- which(index == imax)

    x.imax <- xr[, imax.inds]
    p.imax <- x.imax %*% ginv(x.imax)
    q.imax <- diag(rep(1, n)) - p.imax
    if (step > 1) {
      qs[[step]] <- qs[[step-1]] %*% q.imax
    } else {
      qs[[step]] <- q.imax
    }
    etas[, step] <- p.imax %*% r / w_max

    # Project orthogonally
    r <- q.imax %*% r
    xr[, inactive.inds] <- q.imax %*% xr[, inactive.inds]
    active[step] <- imax

    ## if (is.matrix(Sigma)) {
    ##   added$RSS <- t(r) %*% whitener %*% r
    ##   scale.chisq <- 1
    ## } else {
    ##   added$RSS <- sum(r^2)
    ##   scale.chisq <- Sigma
    ## }

    ## added$RSSdrop <- sum((r.last - r)^2)
    ## added$chisq <- pchisq(added$RSSdrop/scale.chisq, lower.tail=FALSE, df = added$k)
    ## r.last <- r
    ## output <- rbind(output, data.frame(added))
    ## if (verbose) print(added)
  }

  if (verbose) print("Forward stepwise finished, now computing intervals")


  inactive <- labels
  outers <- lapply(labels, function(i) {
    xi <- x[, which(index == i)]
    drop(xi %o% xi)
  })
  aij <- array(NA, dim = c(steps, p, n, n))
  for (step in 1:steps) {
    i <- active[step]
    inactive <- setdiff(inactive, i)
    for (j in inactive) {
      xj <- x[, which(index == j)]
      aij[i, j, ,,] <- outers[i] - outers[j]
    }
  }
  
  if (verbose) print("Formed Aij")

  for (step in 1:steps) {
    v <- active[step]
    eta <- etas[, step]
    as <- bs <- cs <- c()

    if (verbose) print(paste0("Forming intervals for variable: ", v))
    
    inactive <- labels
    for (s in 1:steps) {
      i <- active[s]
      r <- rs[, i]
      inactive <- setdiff(inactive, i)
      for (j in inactive) {
        A <- aij[i, j, ,,]
        etaA <- t(eta) %*% A
        as <- c(as, etaA %*% eta)
        bs <- c(bs, 2 * etaA %*% r)
        cs <- c(cs, t(r) %*% A %*% r)

        # For as > 0, track upper and lower
        # for as < 0, track both parts separately
        # form intersection intelligently
      }
    }
    
  }
  
  ## value <- list(variable=output$imax, p.value=output$p.value, log=output)
  ## class(value) <- "fstep"
  ## attr(value, "n") <- nrow(x)
  ## attr(value, "p") <- ncol(x)
  ## attr(value, "labels") <- labels
  ## attr(value, "index") <- index
  ## attr(value, "steps") <- steps
  ## if (!is.null(attr(x, "varnames"))) {
  ##   attr(value, "varnames") <- attr(x, "varnames")
  ## } else {
  ##   attr(value, "varnames") <- colnames(x)
  ## }
  ## if (!is.matrix(Sigma)) attr(value, "scale") <- Sigma

  invisible(value)
}


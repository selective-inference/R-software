# Functions for computing truncated chi p-value

#' Trigonometry for solving 1-dimensional problem
#'
#' @param num Numerator term
#' @param den Denominator term
#' @param w Weight of chosen group
#' @param tol Numerical tolerance
#' @return Upper and lower limits for truncated chi
#' @references \url{http://arxiv.org/abs/1405.3920}
#' @seealso \code{\link{tchi_cdf}}
tchi_trig_solution = function(num, den, w, tol=1.e-10) {

  norm_num <- sqrt(sum(num^2))
  norm_den <- sqrt(sum(den^2))
  cands <- c()

  # Next two cases: effectively no truncation
  if (norm_den == 0) {
    return(c(0, Inf))
  }

  if ((norm_num / norm_den) < tol) {
    return(c(0, Inf))
  }

  Ctheta <- sum(num * den) / (norm_num * norm_den)
  Ctheta <- min(max(Ctheta, -1), 1)
  theta <- acos(Ctheta)
  # Group size > 1
  if (length(num) > 1) {

    Stheta <- sqrt(1-Ctheta^2)
    Sphi <- norm_den * Stheta / w

    # Again, no truncation:
    if (Sphi > 1) {
      return(c(0, Inf))
    }

    phi1 <- asin(Sphi)
    phi2 <- pi - phi1

  } else { # group size 1
    # asin(Sphi) will return 0, this is captured by V3
    phi1 <- pi
    phi2 <- -pi
    cands = c(cands, norm_num / (w - norm_den * cos(theta)))
  }

  cands <- c(cands, norm_num * cos(phi1) / (w - norm_den * cos(theta-phi1)))
  cands <- c(cands, norm_num * cos(phi2) / (w - norm_den * cos(theta-phi2)))

  if (norm_den < w) {
    # No upper truncation
    return(c(max(cands), Inf))
  }

  #cands <- cands[cands >= 0]
  return(c(max(min(cands), 0), max(cands)))
}


#' Truncated CDF transform
#'
#' @param L Lambda, the maximum of \deqn{\| X_g^T y\|_2}
#' @param lower_bound Truncation lower bound
#' @param upper_bound Truncation upper bound
#' @param cvar Conditional \deqn{\sigma}
#' @param k Rank of maximizing group
#' @return Truncated chi p-value
#' @references \url{http://arxiv.org/abs/1405.3920}
#' @seealso \code{\link{tchi_trig_solution}}
tchi_cdf <- function(L, lower_bound, upper_bound, cvar, k) {

  L2 <- L^2
  lb2 <- lower_bound^2

  if (upper_bound == Inf) {
    first.term <- 1
  } else {
    ub2 <- upper_bound^2
    first.term <- pchisq(ub2/cvar, k, lower.tail=TRUE)
  }

  if (first.term == 1) {
    num <- pchisq(L2/cvar, k, lower.tail=FALSE, log.p=TRUE)
    den <- pchisq(lb2/cvar, k, lower.tail=FALSE, log.p=TRUE)
    value <- exp(num - den)

  } else {
    # Is this numerically unstable?
    num <- first.term - pchisq(L2/cvar, k, lower.tail=TRUE)
    den <- first.term - pchisq(lb2/cvar, k, lower.tail=TRUE)
    value <- num/den
  }

  return(value)
}


#' Monte Carlo estimation of an exact p-value (warning: Sigma and weights not implemented)
#'
#' @param L Lambda, the maximum of \deqn{\| X_g^T y\|_2}
#' @param X Design matrix
#' @param index Indicator vector of group membership
#' @param inactive Labels of inactive groups
#' @param M Number of Monte Carlo samples, default 200
#' @return Monte Carlo estimated exact p-value
MC_pvalue <- function(L, X, index, inactive, M=200) {
  L2 <- L^2
  n = nrow(X)
  inactive.labels <- unique(inactive)
  inactive.inds <- which(index %in% inactive)
  inactive.index <- index[inactive.inds]

  MC_instance <- function() {
    epsilon = rnorm(n)
    U = t(X[, inactive.inds]) %*% epsilon
    terms = lapply(inactive.labels, function(x) sum(U[inactive.index==x]^2))
    return(max(unlist(terms)))
  }

  L2s <- unlist(replicate(M, MC_instance(), simplify = FALSE))
  q = sum(L2s > L2)/M

  return(q)
}


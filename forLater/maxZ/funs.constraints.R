#
# Some utilities for affine constraints
#

#
# compute the square-root and inverse square-root of a non-negative
# definite matrix
#

#' Compute the square-root and inverse square-root of a non-negative definite matrix. 
#' @param S matrix
#' @param rank rank of svd 
#' 
#'                                   
factor_covariance = function(S, rank=NA) {
    if (is.na(rank)) {
        rank = nrow(S)
    }
    svd_X = svd(S, nu=rank, nv=rank)
    sqrt_cov = t(sqrt(svd_X$d[1:rank]) * t(svd_X$u[,1:rank]))
    sqrt_inv = t((1. / sqrt(svd_X$d[1:rank])) * t(svd_X$u[,1:rank]))

    return(list(sqrt_cov=sqrt_cov, sqrt_inv=sqrt_inv))
}

#
# from a constraint, return an equivalent
# constraint and a whitening and inverse
# whitening map
#

# law is Z \sim N(mean_param, covariance) subject to constraints linear_part %*% Z \leq offset

#' Transform non-iid problem into iid problem
#' @param linear_part matrix, linear part of constraints
#' @param offset vector, bias of constraints 
#' @param mean_param vector of unconditional means
#' @param covariance vector of unconditional covariance
#' @return new \code{linear_part} and \code{offset} for 0-mean iid covariance problem, 
#' and functions that map between the two problems.                                   
whiten_constraint = function(linear_part, offset, mean_param, covariance) {

    factor_cov = factor_covariance(as.matrix(covariance))
    sqrt_cov = factor_cov$sqrt_cov
    sqrt_inv = factor_cov$sqrt_inv

    new_A = linear_part %*% sqrt_cov
    new_b = offset - linear_part %*% mean_param

    # rescale rows to have length 1

    scaling = sqrt(apply(new_A^2, 1, sum))
    new_A = new_A / scaling
    new_b = new_b / scaling


    inverse_map = function(Z) {
	# broadcasting here
        # the columns of Z are same length as mean_param
        return(sqrt_cov %*% Z + as.numeric(mean_param))
    }

    forward_map = function(W) {
        return(sqrt_inv %*% (W - mean_param))
    }

    return(list(linear_part=new_A,
                offset=new_b,
                inverse_map=inverse_map,
                forward_map=forward_map))
}

#' Sample from multivariate normal distribution under affine restrictions
#' @description
#' \code{sample_from_constraints} returns a sample from the conditional 
#' multivariate normal Z~ N(mean,covariance) s.t. A*Z <= B 
#'
#' @param linear_part r x d matrix for r restrictions and d dimension of Z
#' @param offset r-dim vector of offsets
#' @param mean_param d-dim mean vector of the unconditional normal
#' @param covariance d x d covariance matrix of unconditional normal
#' @param initial_point d-dim vector that initializes the sampler (must meet restrictions)
#' @param ndraw size of sample 
#' @param burnin samples to throw away before storing
#' @return Z ndraw x d matrix of samples 
#' @export 
#' @examples
#' 
#' truncatedNorm = function(1000, c(0,0,0), identity(3), lower = -1,
#'                  upper = c(2,1,2), start.value = c(0,0,0))
#'
#' constr = thresh2constraints(3, lower = c(1,1,1))
#'
#' samp = sample_from_constraints(linear_part = constr$linear_part,
#'                                    offset= constr$offset,
#'                                    mean_param = c(0,0,0),
#'                                    covariance = diag(3),
#'                                    initial_point = c(1.5,1.5,1.5), 
#'                                    ndraw=100,
#'                                    burnin=2000)
#' 

sample_from_constraints = function(linear_part, 
                                   offset, 
                                   mean_param, 
                                   covariance, 
                                   initial_point, 
                                   ndraw=8000,
                                   burnin=2000) 
{

  whitened_con = whiten_constraint(linear_part,
                                     offset,
                                     mean_param,
                                     covariance)
  white_initial = whitened_con$forward_map(initial_point)


  white_linear = whitened_con$linear_part
  white_offset = whitened_con$offset

	# Inf cannot be used in C code
	# In theory, these rows can be dropped

	rows_to_keep = white_offset < Inf
	white_linear = white_linear[rows_to_keep,,drop=FALSE]
	white_offset = white_offset[rows_to_keep]

        nstate = length(white_initial)
	if (sum(rows_to_keep) > 0) {
	    if (ncol(white_linear) > 1) {
                nconstraint = nrow(white_linear)

                directions = rbind(diag(rep(1, nstate)),
                                   matrix(rnorm(nstate^2), nstate, nstate))

                # normalize rows to have length 1

                scaling = apply(directions, 1, function(x) {  return(sqrt(sum(x^2))) })
                directions = directions / scaling
                ndirection = nrow(directions)

                alphas = directions %*% t(white_linear)
                U = white_linear %*% white_initial - white_offset
                Z_sample = matrix(rep(0, nstate * ndraw), nstate, ndraw)

                result = .C("sample_truncnorm_white",
                            as.numeric(white_initial),
                            as.numeric(U),
                            as.numeric(t(directions)),
                            as.numeric(t(alphas)),
                            output=Z_sample,
                            as.integer(nconstraint),
                            as.integer(ndirection),
                            as.integer(nstate),
                            as.integer(burnin),
                            as.integer(ndraw),
                            package="selectiveInference")
                Z_sample = result$output
           } else { # the distribution is univariate
                    # we can just work out upper and lower limits

                white_linear = as.numeric(white_linear)
                pos = (white_linear * white_offset) >= 0
                neg = (white_linear * white_offset) <= 0
		if (sum(pos) > 0) {
                    U = min((white_offset / white_linear)[pos])
                } else {
		    U = Inf
                }
		if (sum(neg) < 0) {
                    L = max((white_offset / white_linear)[neg])
                } else {
                    L = -Inf
                }
                Z_sample = matrix(qnorm((pnorm(U) - pnorm(L)) * runif(ndraw) + pnorm(L)), 1, ndraw)
           }
       } else {
           Z_sample = matrix(rnorm(nstate * ndraw), nstate, ndraw)
  }


  Z = t(whitened_con$inverse_map(Z_sample))
  return(Z)
}

#' Translate between coordinate thresholds and affine constraints
#' @description
#' \code{thresh2constraints} translates lower and upper constraints 
#' on coordinates into linear and offset constraints (A*Z <= B).
#' lower and upper can have -Inf or Inf coordinates.
#' @param d dimension of vector
#' @param lower 1 or d-dim lower constraints 
#' @param upper 1 or d-dim upper constraints 
#' @export
thresh2constraints = function(d, lower = rep(-Inf, d), upper = rep(Inf,d)){
  stopifnot(is.element(length(lower),c(1,d)))
  stopifnot(is.element(length(upper),c(1,d)))

  if (length(lower) == 1){
    lower = rep(lower, d)
  }
  if (length(upper) == 1){
    upper = rep(upper, d)
  }


  linear_part = matrix(ncol = d, nrow = 0)
  offset = numeric(0)
  lower_constraints = which(lower > -Inf)
  for (l in lower_constraints){
    new_vec = rep(0,d)
    new_vec[l] = -1
    linear_part = rbind(linear_part, new_vec)
    offset = c(offset, -lower[l])
  }
  upper_constraints = which(upper < Inf)
  for (u in upper_constraints){
    new_vec = rep(0,d)
    new_vec[u] = 1
    linear_part = rbind(linear_part, new_vec)
    offset = c(offset, upper[u])
  }

  constraints = list(linear_part = linear_part, offset = offset)
  return(constraints)
}




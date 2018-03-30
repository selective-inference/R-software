#
# Some utilities for affine constraints
#

# 
# compute the square-root and inverse square-root of a non-negative
# definite matrix
#

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

whiten_constraint = function(linear_part, offset, mean_param, covariance) {

    factor_cov = factor_covariance(covariance)
    sqrt_cov = factor_cov$sqrt_cov
    sqrt_inv = factor_cov$sqrt_inv

    new_A = linear_part %*% sqrt_cov
    new_b = offset - linear_part %*% mean_param

    # rescale rows to have length 1

    scaling = sqrt(apply(new_A^2, 1, sum))
    new_A = new_A / scaling
    new_b = new_b / scaling

    # TODO: check these functions will behave when Z is a matrix.

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

#
# sample from the law
#
# Z \sim N(mean_param, covariance) subject to constraints linear_part %*% Z \leq offset

sample_from_constraints = function(linear_part,
                                   offset,
                                   mean_param,
                                   covariance,
                                   initial_point, # point must be feasible for constraints
                                   ndraw=8000,
                                   burnin=2000,
                                   accept_reject_params=NA) #TODO: implement accept reject check
{

    whitened_con = whiten_constraint(linear_part,
                                     offset,
				     mean_param,
                                     covariance)
    white_initial = whitened_con$forward_map(initial_point)

#     # try 100 draws of accept reject
#     # if we get more than 50 good draws, then just return a smaller sample
#     # of size (burnin+ndraw)/5

#     if accept_reject_params: 
#         use_hit_and_run = False
#         num_trial, min_accept, num_draw = accept_reject_params

#         def _accept_reject(sample_size, linear_part, offset):
#             Z_sample = np.random.standard_normal((100, linear_part.shape[1]))
#             constraint_satisfied = (np.dot(Z_sample, linear_part.T) - 
#                                     offset[None,:]).max(1) < 0
#             return Z_sample[constraint_satisfied]

#         Z_sample = _accept_reject(100, 
#                                   white_con.linear_part,
#                                   white_con.offset)

#         if Z_sample.shape[0] >= min_accept:
#             while True:
#                 Z_sample = np.vstack([Z_sample,
#                                       _accept_reject(num_draw / 5,
#                                                      white_con.linear_part,
#                                                      white_con.offset)])
#                 if Z_sample.shape[0] > num_draw:
#                     break
#             white_samples = Z_sample
#         else:
#             use_hit_and_run = True
#     else:
#         use_hit_and_run = True

    use_hit_and_run = TRUE

    if (use_hit_and_run) {

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
    }

    Z = t(whitened_con$inverse_map(Z_sample))
    return(Z)
}


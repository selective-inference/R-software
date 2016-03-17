#
# Some utilities for affine constraints
#

# 
# compute the square-root and inverse square-root of a non-negative
# definite matrix
#

factor_covariance = function(S, rank=NA) {
    if (!is.na(rank)) {
        rank = nrow(S)
    }
    svd_X = svd(S, nu = rank, nv=rank)
    sqrt_cov = t(sqrt(svd_X$d[1:rank]) * t(svd_X$u[,1:rank]))
    sqrt_inv = t((1. / sqrt(svd_X$d[1:rank])) * t(svd_X$u[,1:rank]))

    return(list(sqrt_cov=sqrt_cov, sqrt_inv=sqrt_inv))
}

#
# from a constraint, return an equivalent
# constraint and a whitening and inverse
# whitening map
#

# law is Z \sim N(mean, covariance) subject to constraints linear_part %*% Z \leq offset

whiten_constraint = function(linear_part, offset, mean, covariance) {

    factor_cov = factor_covariance(covariance)
    sqrt_cov = factor_cov$sqrt_cov
    sqrt_inv = factor_cov$sqrt_inv

    new_A = A %*% sqrt_cov
    new_b = offset - linear_part %*% mean

    # rescale rows to have length 1

    scaling = sqrt(apply(new_A^2, sum, 1))
    new_A = new_A / scaling
    new_b = new_b / scaling

    # TODO: check these functions will behave when Z is a matrix.

    inverse_map = function(Z) {
        return(sqrt_cov %*% Z + mu)
    }

    forward_map = function(W) {
        return(sqrt_inv %*% (W - mu))
    }    

    return(list(linear_part=new_A,
                offset=new_b,
                inverse_map=inverse_map,
                forward_map=forward_map))
}

#
# sample from the law
#
# Z \sim N(mean, covariance) subject to constraints linear_part %*% Z \leq offset

sample_from_constraints = function(linear_part,
                                   offset,
                                   covariance,
                                   mean,
                                   initial_point, # point must be feasible for constraints
                                   ndraw=8000,
                                   burnin=2000,
                                   accept_reject_params=NA) #TODO: implement accept reject check
{

    whitened_con = whiten_constraint(linear_part,
                                     offset,
                                     covariance,
                                     mean)
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

        nstate = length(white_initial)
        directions = rbind(diag(rep(1, ndim)),
                           matrix(rnorm(ndim^2), nstate, nstate))
        directions = apply(directions, function(x) { x/sqrt(sum(x^2)) }, 1) # normalize rows to have length 1
        alphas = white_linear %*% directions
        U = white_linear %*% white_initial - white_offset
        Z_sample = matrix(rep(NA, ndraw*ndim), nstate, ndraw)

        .C("sample_truncnorm_white",
           white_initial,
           U,
           t(directions),
           t(alphas),
           Z_sample,
           nrow(white_linear),
           nstate,
           burnin,
           ndraw)
    
    }
    Z = t(whitened_con$inverse_map(Z_sample))
    return(Z)
}


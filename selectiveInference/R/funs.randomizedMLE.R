# This function computes the approximate UMVU as described
# in Lemma 4 of https://arxiv.org/abs/1703.06176

# - \Lambda^*_g is the quadratic form with randomizer_precision
# - P is opt_transform$linear_term
# - D is target_transform$linear_term
# - q is target_transform$offset_term + opt_transform$offset_term

solve_UMVU = function(target_transform,
                      opt_transform,
                      target_observed,
                      feasible_point,
                      target_cov,
                      randomizer_precision,
                      nstep=30,
                      tol=1.e-8) {

    D = target_transform$linear_term; data_offset = target_transform$offset_term
    P = opt_transform$linear_term; opt_offset = opt_transform$offset_term
    q = data_offset + opt_offset

    nopt = ncol(P)
    ntarget = ncol(D)
    ntotal = nopt + ntarget

    # setup joint implied covariance matrix

    target_precision = solve(target_cov)
    implied_precision = matrix(0, ntarget + nopt, ntarget + nopt)

    implied_precision[1:ntarget,1:ntarget] = t(D) %*% randomizer_precision %*% D + target_precision
    implied_precision[1:ntarget,(ntarget+1):ntotal] = t(D) %*% randomizer_precision %*% P
    implied_precision[(ntarget+1):ntotal,1:ntarget] = t(P) %*% randomizer_precision %*% D
    implied_precision[(ntarget+1):ntotal,(ntarget+1):ntotal] = t(P) %*% randomizer_precision %*% P
    implied_cov = solve(implied_precision)

    implied_opt = implied_cov[(ntarget+1):ntotal,(ntarget+1):ntotal,drop=FALSE]
    implied_target = implied_cov[1:ntarget,1:ntarget,drop=FALSE]
    implied_cross = implied_cov[1:ntarget,(ntarget+1):ntotal,drop=FALSE]

    L = implied_cross %*% solve(implied_opt)
    M_1 = solve(implied_precision[1:ntarget,1:ntarget,drop=FALSE]) %*% target_precision
    M_2 = -solve(implied_precision[1:ntarget,1:ntarget,drop=FALSE]) %*% t(D) %*% randomizer_precision

    conditioned_value = q

    linear_term = t(implied_cross) %*% solve(implied_target)
    offset_term = -t(P) %*% randomizer_precision %*% conditioned_value
    natparam_transform = list(linear_term=linear_term, offset_term=offset_term)
    conditional_natural_parameter = linear_term %*% target_observed + offset_term

    conditional_precision = implied_precision[(ntarget+1):ntotal,(ntarget+1):ntotal,drop=FALSE]

    M_1_inv = solve(M_1)
    offset_term = - M_1_inv %*% M_2 %*% conditioned_value
    mle_transform = list(target_lin=M_1_inv, soln_lin=-M_1_inv.dot(L), offset=offset_term)

    mle_map = function(target_observed, feasible_point=rep(1, length(target_observed))) {
        param_lin = natparam_transform$linear_term
	param_offset = natparam_transform$offset_term
        mle_target_lin = mle_transfrom$target_lin
	mle_soln_lin = mle_transform$soln_lin
	mle_offset = mle_transform$offset

        result = solve_barrier_(param_lin %*% target_observed + param_offset,
                                conditional_precision,
                                feasible_point,
                                nstep,
                                tol,
                                mean(diag(conditional_precision)))

        return(list(soln=mle_target_lin %*% target_observed + mle_soln_lin %*% soln + mle_offset, 
                    value=value))
    }
    sel_MLE = mle_map(target_observed)
    return(list(soln=sel_MLE$soln, value=sel_MLE$value, map=mle_map))
}

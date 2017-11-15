library(selectiveInference)

simple_problem = function(target_observed=-3., threshold=2., randomization_scale=1.) {
    # Simple problem: randomizaiton of sd 1 and thresholded at 2 (default args)

    target_observed = as.numeric(target_observed)
    target_transform = list(linear_term=-diag(1), offset_term=c(0.))
    opt_transform = list(linear_term=diag(1), offset_term=c(threshold*1.))
    feasible_point = c(1.)
    randomizer_precision = diag(1. / randomization_scale^2)
    target_cov = diag(1)

    print(target_transform)
    print(opt_transform)
    print(feasible_point)
    print(target_cov)
    return(selectiveInference:::solve_UMVU(target_transform,
                                           opt_transform,
                                           target_observed,
                                           feasible_point,
                                           target_cov,
                                           randomizer_precision))
}

simple_problem()

library(selectiveInference)

smoke_test = function() {

    n = 100; p = 50
    X = matrix(rnorm(n * p), n, p)
    y = rnorm(n)
    lam = 20 / sqrt(n)
    noise_scale = 0.01 * sqrt(n)
    ridge_term = .1 / sqrt(n)
    selectiveInference:::randomizedLASSO(X, y, lam, noise_scale, ridge_term)
}
A = smoke_test()

density_test = function() {

    random_lasso = smoke_test()
    p = nrow(random_lasso$internal_transform$linear_term)
    internal_state = matrix(rnorm(p * 20), p, 20)
    optimization_state = matrix(rnorm(p * 20), p, 20)
    offset = rnorm(p)

    selectiveInference:::log_density_gaussian_(10.,
                                               random_lasso$internal_transform$linear_term,
                                               internal_state,
                                               random_lasso$optimization_transform$linear_term,
					       optimization_state,
                                               offset)

    selectiveInference:::log_density_gaussian_conditional_(10.,
                                                           random_lasso$optimization_transform$linear_term,
                    				           optimization_state,
                                                           offset)
}

density_test()

library(selectiveInference)

smoke_test = function() {
    n = 100; p = 50
    X = matrix(rnorm(n * p), n, p)
    y = rnorm(n)
    lam = 20 / sqrt(n)
    noise_scale = 0.01 * sqrt(n)
    ridge_term = .1 / sqrt(n)
    selectiveInference:::randomizedLasso(X, y, lam, noise_scale, ridge_term)
}

A = smoke_test()

sampler_test = function() {

    n = 100; p = 50
    X = matrix(rnorm(n * p), n, p)
    y = rnorm(n)
    lam = 20 / sqrt(n)
    noise_scale = 0.01 * sqrt(n)
    ridge_term = .1 / sqrt(n)
    obj = selectiveInference:::randomizedLasso(X, y, lam, noise_scale, ridge_term)
    S = selectiveInference:::sample_opt_variables(obj, jump_scale=rep(1/sqrt(n), p), nsample=10000)
    return(S$samples[2001:10000,])
}
B = sampler_test()

gaussian_density_test = function() {

    noise_scale = 10.
    random_lasso = smoke_test()
    p = nrow(random_lasso$internal_transform$linear_term)
    internal_state = matrix(rnorm(p * 20), p, 20)
    optimization_state = matrix(rnorm(p * 20), p, 20)
    offset = rnorm(p)

    V1 = selectiveInference:::log_density_gaussian_(noise_scale,
                                                    random_lasso$internal_transform$linear_term,
                                                    internal_state,
                                                    random_lasso$optimization_transform$linear_term,
                                                    optimization_state,
                                                    offset)
    A1 = random_lasso$internal_transform$linear_term
    A2 = random_lasso$optimization_transform$linear_term
    arg = A1 %*% internal_state + A2 %*% optimization_state + offset
    V2 = -apply(arg^2, 2, sum) / (2 * noise_scale^2)
    print(sqrt(sum((V1-V2)^2) / sum(V1^2)))

    U1 = selectiveInference:::log_density_gaussian_conditional_(noise_scale,
                                                                random_lasso$optimization_transform$linear_term,
                                                                optimization_state,
                                                                offset)
    arg = A2 %*% optimization_state + offset
    U2 = -apply(arg^2, 2, sum) / (2 * noise_scale^2)
    print(sqrt(sum((U1-U2)^2) / sum(U1^2)))

    # test that a single column matrix works -- numeric should not

    print(selectiveInference:::log_density_gaussian_conditional_(noise_scale,
                                                                 random_lasso$optimization_transform$linear_term,
                                                                 optimization_state[,1,drop=FALSE],
                                                                 offset))
    print(selectiveInference:::log_density_gaussian_(noise_scale,
                                                     random_lasso$internal_transform$linear_term,
                                                     internal_state[,1,drop=FALSE],
                                                     random_lasso$optimization_transform$linear_term,
                                                     optimization_state[,1,drop=FALSE],
                                                     offset))

}

gaussian_density_test()

laplace_density_test = function() {

    noise_scale = 10.
    random_lasso = smoke_test()
    p = nrow(random_lasso$internal_transform$linear_term)
    internal_state = matrix(rnorm(p * 20), p, 20)
    optimization_state = matrix(rnorm(p * 20), p, 20)
    offset = rnorm(p)

    V1 = selectiveInference:::log_density_laplace_(noise_scale,
                                                    random_lasso$internal_transform$linear_term,
                                                    internal_state,
                                                    random_lasso$optimization_transform$linear_term,
                                                    optimization_state,
                                                    offset)
    A1 = random_lasso$internal_transform$linear_term
    A2 = random_lasso$optimization_transform$linear_term
    arg = A1 %*% internal_state + A2 %*% optimization_state + offset
    V2 = -apply(abs(arg), 2, sum) / noise_scale
    print(sqrt(sum((V1-V2)^2) / sum(V1^2)))

    U1 = selectiveInference:::log_density_laplace_conditional_(noise_scale,
                                                                random_lasso$optimization_transform$linear_term,
                                                                optimization_state,
                                                                offset)
    arg = A2 %*% optimization_state + offset
    U2 = -apply(abs(arg), 2, sum) / noise_scale
    print(sqrt(sum((U1-U2)^2) / sum(U1^2)))

    # test that a single column matrix works -- numeric should not

    print(selectiveInference:::log_density_laplace_conditional_(noise_scale,
                                                                 random_lasso$optimization_transform$linear_term,
                                                                 optimization_state[,1,drop=FALSE],
                                                                 offset))
    print(selectiveInference:::log_density_laplace_(noise_scale,
                                                     random_lasso$internal_transform$linear_term,
                                                     internal_state[,1,drop=FALSE],
                                                     random_lasso$optimization_transform$linear_term,
                                                     optimization_state[,1,drop=FALSE],
                                                     offset))

}

laplace_density_test()

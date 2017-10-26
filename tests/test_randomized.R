library(selectiveInference)

test = function() {

    n = 100; p = 50
    X = matrix(rnorm(n * p), n, p)
    y = rnorm(n)
    lam = 20 / sqrt(n)
    noise_scale = 0.01 * sqrt(n)
    ridge_term = .1 / sqrt(n)
    fit_randomized_lasso(X, y, lam, noise_scale, ridge_term)
}

print(test())

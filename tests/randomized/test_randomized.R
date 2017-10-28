library(selectiveInference)

test = function() {

    n = 100; p = 50
    X = matrix(rnorm(n * p), n, p)
    y = rnorm(n)
    lam = 20 / sqrt(n)
    noise_scale = 0.01 * sqrt(n)
    ridge_term = .1 / sqrt(n)
    selectiveInference:::randomizedLASSO(X, y, lam, noise_scale, ridge_term)
}

A=test()
#print(test())

library(selectiveInference)
### Test

n = 100; p = 50

X = matrix(rnorm(n * p), n, p)
Y = rnorm(n)
lam = 2

soln1 = selectiveInference:::fit_randomized_lasso(X, Y, lam, 0, 0)$soln
G = glmnet(X, Y, intercept=FALSE, standardize=FALSE)
soln2 = coef(G, s=1/n, exact=TRUE, x=X, y=Y)[-1]

print(soln1)
print(soln2)
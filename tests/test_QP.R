library(selectiveInference)
### Test
n = 80; p = 50


X = matrix(rnorm(n * p), n, p)
Y = rnorm(n)
lam = 2

soln1 = selectiveInference:::fit_randomized_lasso(X, Y, lam, 0, 0)$soln
G = glmnet(X, Y, intercept=FALSE, standardize=FALSE)
soln2 = coef(G, s=lam/n, exact=TRUE, x=X, y=Y)[-1]

print(soln1)
print(soln2)
plot(soln1, soln2)
print(summary(lm(soln1 ~ soln2)))


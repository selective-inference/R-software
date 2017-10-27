library(selectiveInference)
source('oldcode.R')

n = 500; p = 50

X = matrix(rnorm(n * p), n, p)
S = t(X) %*% X / n

mu = 7.791408e-02

A1 = debiasingMatrix(S, FALSE, n, 1:5, mu=mu, max_iter=1000)
A2 = debiasingMatrix(S / n, FALSE, n, 1:5, mu=mu, max_iter=1000)

B1 = debiasingMatrix(X, TRUE, n, 1:5, mu=mu, max_iter=1000)
B2 = debiasingMatrix(X / sqrt(n), TRUE, n, 1:5, mu=mu, max_iter=1000)

C1 = InverseLinfty(S, n, mu=mu, maxiter=1000)[1:5,]
C2 = InverseLinfty(S / n, n, mu=mu, maxiter=1000)[1:5,]

par(mfrow=c(2,3))
plot(A1[1,], C1[1,])
plot(A1[1,], B1[1,])
plot(B1[1,], C1[1,])

plot(A1[1,], A2[1,])
plot(B1[1,], B2[1,])
plot(C1[1,], C2[1,])

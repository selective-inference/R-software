library(lars)
library(polypath)

set.seed(0)
n = 20
p = 25
y = rnorm(n)
X = matrix(rnorm(n*p),n,p)

# LAR
out1 = lars(X,y,type="lar",normalize=FALSE,intercept=FALSE)
out2 = lar(y,X,verbose=TRUE)

# Some checks against the lars package
max(abs(out1$lambda-out2$lambda))
max(abs(t(out1$beta[-nrow(out1$beta),])-out2$beta))

# Some checks for the Gamma matrix
lam2 = out2$Gamma[out2$nk,] %*% y
max(abs(out2$lambda-lam2))
min(out2$Gamma %*% y)
nrow(out2$Gamma)

# Lasso
out3 = lars(X,y,type="lasso",normalize=FALSE,intercept=FALSE)
out4 = lasso(y,X,verbose=TRUE)

# Some checks against the lars package
max(abs(out3$lambda-out4$lambda))
max(abs(t(out3$beta[-nrow(out3$beta),])-out4$beta))

# Some checks for the Gamma matrix
lam2 = out4$Gamma[out4$nk,] %*% y
max(abs(out4$lambda-lam2))
min(out4$Gamma %*% y)
nrow(out4$Gamma)


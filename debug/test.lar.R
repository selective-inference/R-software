library(selectiveInference)
library(lars)

set.seed(0)
n = 25
p = 50
s = 3
size = 10

sigma = 1
x = matrix(rnorm(n*p),n,p)
b = c(sample(c(-1,1),s,replace=T)*rep(size,s),rep(0,p-s))
mu = x%*%b
y = mu + sigma*rnorm(n)

obj = lar(x,y,verb=T,intercept=T,norm=T)
obj2 = lars(x,y,intercept=T,norm=T,type="lar")

# Checks
max(abs(obj$lambda - obj$Gamma[obj$nk,] %*% obj$y))

max(abs(obj$lambda - obj2$lambda))

max(abs(coef.lar(obj,s=(obj$lambda[4]+obj$lambda[5])/2,mode="lam")-
        lars::predict.lars(obj2,s=(obj$lambda[4]+obj$lambda[5])/2,type="coef",mode="lam")$coef))

max(abs(predict.lar(obj,s=4.5,mode="step")-
        lars::predict.lars(obj2,s=4.5,newx=x,mode="step")$fit))

# Sequential inference
out = larInf(obj,sigma=sigma)

if (false) {
plot(out$pv,ylim=c(0,1))
points(out$pv.spacing,col=2)
points(out$pv.asymp,col=3)
points(out$pv.covtest,col=4)
legend("topleft",
       legend=c("Exact","Spacing","Asymp","Cov test"),
       pch=21,col=1:4)
}

a = lm(y~x[,1:s]+0)
stdcoefs = a$coef / sqrt(diag(vcov(a)))

round(cbind(out$ci, out$vmat %*% mu * out$sign, c(stdcoefs,rep(0,nrow(out$ci)-s))),3)
sum(out$ci[,1]>out$ci[,2])

# Fixed step inference
k = 8
out = larInf(obj,sigma=sigma,k=k,type="all")

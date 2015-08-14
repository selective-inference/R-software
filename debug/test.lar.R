#library(selectiveInference)
library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")
library(lars)

set.seed(0)
n = 50
p = 100
s = 3
size = 3

sigma = .2
x = matrix(rnorm(n*p),n,p)
b = c(sample(c(-1,1),s,replace=T)*rep(size,s),rep(0,p-s))
mu = x%*%b
y = mu + sigma*rnorm(n)

obj = lar(x,y,verb=T,intercept=T,norm=T)
obj2 = lars(x,y,intercept=T,norm=T,type="lar")

plot(obj)

# Checks
max(abs(obj$lambda - obj$Gamma[obj$nk,] %*% obj$y))
max(abs(obj$lambda - obj2$lambda))
max(abs(coef(obj,s=(obj$lambda[4]+obj$lambda[5])/2,mode="lam")-
        lars::predict.lars(obj2,s=(obj$lambda[4]+obj$lambda[5])/2,type="coef",mode="lam")$coef))
max(abs(predict(obj,s=4.5,mode="step")-
        lars::predict.lars(obj2,s=4.5,newx=x,mode="step")$fit))

# Sequential inference
out = larInf(obj,sigma=sigma,k=20,verbose=TRUE)
out
sum(out$ci[,1]>out$ci[,2])
plot(out$pv,ylim=c(0,1))

# AIC inference
k = 20
out2 = larInf(obj,sigma=sigma,k=k,type="aic")
out2

# Fixed step inference
k = out2$khat
out3 = larInf(obj,sigma=sigma,k=k,type="all")
print(out3,tail=TRUE)

# Least squares inference
out.ls = lm(y~x[,obj$action[1:k]])
summary(out.ls)

# Don't lose much, in terms of conditioning on AIC event,
# But the p-values don't look great here ...
# We don't get significance for variables 1 and 2. When you
# take k=5 steps, we do see significance

# Fixed step inference
larInf(obj,sigma=sigma,k=4,type="all")

### It seems like the presence of other variables in the model,
### under conditioning, messes with the p-values
### In other words, correlation between variables affects
### variable significance MORE with conditioning that without

#################
#################
# Another random seed, a little more favorable results

set.seed(1)
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

# Sequential inference
out = larInf(obj,sigma=sigma)
out

# AIC  inference
k = 15
out2 = larInf(obj,sigma=sigma,k=k,type="aic")
out2

# Fixed step inference
k = out2$khat
out3 = larInf(obj,sigma=sigma,k=k,type="all")
out3

# Least squares inference
out.ls = lm(y~x[,obj$action[1:k]])
summary(out.ls)

# Explore fixed step inferences
larInf(obj,sigma=sigma,k=3,type="all")
larInf(obj,sigma=sigma,k=4,type="all")
larInf(obj,sigma=sigma,k=5,type="all")
larInf(obj,sigma=sigma,k=6,type="all")
larInf(obj,sigma=sigma,k=7,type="all")
larInf(obj,sigma=sigma,k=8,type="all")
larInf(obj,sigma=sigma,k=9,type="all")
larInf(obj,sigma=sigma,k=10,type="all")

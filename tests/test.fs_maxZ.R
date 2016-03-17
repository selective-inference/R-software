library(selectiveInference)
options(error=dump.frames)

set.seed(0)
n = 100
p = 120
s = 3
size = 5

sigma = 1.5
x = matrix(rnorm(n*p),n,p)
#x = scale(x,T,F)/sqrt(n-1)

b = c(sample(c(-1,1),s,replace=T)*rep(size,s),rep(0,p-s))
mu = x%*%b
y = mu + sigma*rnorm(n)

obj = fs(x,y,verb=T,intercept=T,norm=T, maxsteps=20)


# Sequential inference
out = fsInf_maxZ(obj,sigma=sigma,k=20, ndraw=5000, burnin=1000)
print(out)

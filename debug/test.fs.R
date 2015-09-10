library(selectiveInference)
#library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

set.seed(0)
n = 100
p = 100
s = 3
size = 5

sigma = 1
x = matrix(rnorm(n*p),n,p)
#x = scale(x,T,F)/sqrt(n-1)

b = c(sample(c(-1,1),s,replace=T)*rep(size,s),rep(0,p-s))
mu = x%*%b
y = mu + sigma*rnorm(n)

obj = fs(x,y,verb=T,intercept=T,norm=T)

# NOTE this does not line up with lars' stepwise function,
# but that's OK, because they used different update rules

# Sequential inference
out = fsInf(obj,sigma=sigma,k=20)
out
sum(out$ci[,1]>out$ci[,2])
plot(out$pv,ylim=c(0,1))

# AIC inference
k = 20
out2 = fsInf(obj,sigma=sigma,k=k,type="aic")
out2

# Fixed step inference
k = out2$khat
out3 = fsInf(obj,sigma=sigma,k=k,type="all")
out3

# Least squares inference
X = x[,obj$action[1:k]]
out.ls = lm(y~X+0)
summary(out.ls)

# Don't lose much, in terms of conditioning on AIC event,
# The p-values look good here! 

#################
#################
# Another random seed

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



#check coverage
set.seed(32)

n=50
p=10
sigma=2

x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,T)/sqrt(n-1)    #try with and without standardization

beta=c(5,4,3,2,1,rep(0,p-5))

nsim=100
seeds=sample(1:9999,size=nsim)
pv=rep(NA,nsim)
ci=matrix(NA,nsim,2)
btrue=rep(NA,nsim)
  mu=x%*%beta
for(ii in 1:nsim){
    cat(ii)
    set.seed(seeds[ii])
  
   y=mu+sigma*rnorm(n)
    y=y-mean(y)  
   fsfit=fs(x,y,norm=T)
  
     junk= fsInf(fsfit,sigma=sigma)
    pv[ii]=junk$pv[1]
    oo=junk$var[1]
     btrue[ii]=lsfit(x[,oo],mu)$coef[2]
     ci[ii,]=junk$ci[1,]
}

sum(ci[,1]> btrue)
sum(ci[,2]< btrue)





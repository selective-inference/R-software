library(selectiveInference)

options(error=dump.frames)


library(lars)

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
obj2 = lars(x,y,type="step",intercept=T,norm=T)

max(abs(obj$action-unlist(obj2$action)))
# These don't always match ... what is the lars function doing?

# Checks
max(abs(obj$action-unlist(obj2$action)))
max(abs(coef(obj,s=4.5,mode="step")-
        lars::predict.lars(obj2,s=4.5,type="coef",mode="step")$coef))
max(abs(predict(obj,s=4.5,mode="step")-
        lars::predict.lars(obj2,s=4.5,newx=x,mode="step")$fit))

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
beta=rep(0,p)
nsim=500
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
plot((1:nsim)/nsim,sort(pv))
    abline(0,1)
    
    
sum(ci[,1]> btrue)
sum(ci[,2]< btrue)



##diabetes example
    x=read.table("/Users/tibs/dropbox/PAPERS/FourOfUs/data64.txt")
x=as.matrix(x)
x=scale(x,T,F)
#x=scale(x,T,T)
n=length(y)

    nams=scan("/Users/tibs/dropbox/PAPERS/FourOfUs/data64.names",what="")
y=scan("/Users/tibs/dropbox/PAPERS/FourOfUs/diab.y")
y=y-mean(y)

obj = fs(x,y,verb=T,intercept=T,norm=T)

# Sequential inference

 sigma= estimateSigma(x,y)$sigmahat
out = fsInf(obj,sigma=sigma,k=20)
out


# AIC inference

out2 = fsInf(obj,sigma=sigma,type="aic")
out2

# Fixed step inference
k = out2$khat
out3 = fsInf(obj,sigma=sigma,k=k,type="all")
out3
    out4 = fsInf(obj,sigma=sigma,k=k,type="all",bits=200)
    
##plot

    library(selectiveInference)

options(error=dump.frames)


     set.seed(33)
     n = 50
     p = 10
     sigma = 1
     x = matrix(rnorm(n*p),n,p)
     beta = c(3,2,rep(0,p-2))
     y = x%*%beta + sigma*rnorm(n)
     
     # run forward stepwise, plot results
     fsfit = fs(x,y)
     plot(fsfit)

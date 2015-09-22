library(selectiveInference)

library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")

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
max(abs(obj$action-unlist(obj2$action))
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
    
###
    
myfs=function(x,y,nsteps=min(nrow(x),ncol(x)),mode=c("ip","ip2","cor","cor2","rss")){
p=ncol(x)
ip=rep(NA,p)
# fs by minimizing scaled ip
# first center x and y
x=scale(x,T,F)
y=y-mean(y)
rss0=sum( (y-mean(y))^2)
pred=s=scor=bhat=rep(NA,nsteps)
  if(mode=="ip"| mode=="ip2")   ip=t(x)%*%y/sqrt(diag(t(x)%*%x))
if(mode=="cor" | mode=="cor2") ip=abs(cor(x,y))

if(mode=="rss") {for(j in 1:p){
                a=lsfit(x[,j],y)
                ip[j]=rss0-sum(a$res^2)
            }}
  pred[1]=which.max(abs(ip))
 s[1]=sign(sum(x[,pred[1]]*y))
scor[1]=ip[pred[1]]
bhat[1]=ip[pred[1]]/sqrt(sum(x[,pred[1]]^2))

r=lsfit(x[,pred[1]],y)$res
for(j in 2:nsteps){
    cat(j)
  mod=pred[1:(j-1)]
  r= lsfit(x[,mod],r)$res
  xr= lsfit(x[,mod],x)$res
 if(mode=="ip")  ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
     if(mode=="ip2")  ip=t(x)%*%r/sqrt(diag(t(x)%*%x))
  if(mode=="cor") ip=abs(cor(xr,r))
      if(mode=="cor2") ip=abs(cor(x,r))
  if(mode=="rss"){
        for(k in (1:p)[-mod]){
                a=lsfit(x[,c(mod,k)],r)
                ip[k]=rss0-sum(a$res^2)
            }}
           
 ip[mod]=0
  pred[j]=which.max(abs(ip))
  scor[j]=ip[pred[j]]
  s[j]=sign(sum(xr[,pred[j]]*r))
 bhat[j]=ip[pred[j]]/sqrt(sum(xr[,pred[j]]^2))
}
return(list(pred=pred,s=s,scor=scor,bhat=bhat))
}


    
set.seed(144)
n <- 40
p <- 80
maxsteps <- 10
sparsity <- 5
snr <- 3
x=matrix(rnorm(n*p),n,p)
#    x=scale(x,T,T)
# Compare to step function in R
index <- 1:ncol(x)
y <- rnorm(n)
beta <- rep(0, p)
beta[which(index %in% 1:sparsity)] <- snr
y <- y + x %*% beta
df <- data.frame(y = y, x = x)
fsfit <- step(lm(y ~ 1, df), direction="forward", scope = formula(lm(y~., df)), steps = 20)

junk2=lars(x,y,type="step")
    a=fs(x,y)
    a2=myfs(x,y,mode="ip",nsteps=20)
    a22=myfs(x,y,mode="ip2",nsteps=20)
    a3=myfs(x,y,mode="cor",nsteps=20)
    a33=myfs(x,y,mode="cor2",nsteps=20)
    a4=myfs(x,y,mode="rss",nsteps=20)
out=cbind(unlist(junk2$act)[1:20],names(fsfit$coefficients)[-1],a$act[1:20],a2$pred,a22$pred,a3$pred,a33$pred,a4$pred)
colnames(out)=c("lars","step","myfs","myfs/ip","myfs/ip2","myfs/cor","myfs/cor2","myfs/rss")

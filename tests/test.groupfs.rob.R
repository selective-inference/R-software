#library(lars)
#library(intervals)
#source("../selectiveInference/R/funs.common.R")
#source("../selectiveInference/R/funs.groupfs.R")
source("../selectiveInference/R/funs.quadratic.R")
##source("../selectiveInference/R/funs.fs.R")
#source("../selectiveInference/R/funs.lar.R")

library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")



myfs=function(x,y,nsteps=min(nrow(x),ncol(x)),mode=c("ip","cor")){
p=ncol(x)
# fs by minimizing scaled ip
# first center x and y
x=scale(x,T,F)
y=y-mean(y)
pred=s=scor=bhat=rep(NA,nsteps)
  if(mode=="ip")   ip=t(x)%*%y/sqrt(diag(t(x)%*%x))
if(mode=="cor") ip=abs(cor(x,y))
  pred[1]=which.max(abs(ip))
 s[1]=sign(sum(x[,pred[1]]*y))
scor[1]=ip[pred[1]]
bhat[1]=ip[pred[1]]/sqrt(sum(x[,pred[1]]^2))

r=lsfit(x[,pred[1]],y)$res
for(j in 2:nsteps){
  mod=pred[1:(j-1)]
  r= lsfit(x[,mod],r)$res
  xr= lsfit(x[,mod],x)$res
 if(mode=="ip")  ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
  if(mode=="cor") ip=abs(cor(xr,r))
 ip[mod]=0
  pred[j]=which.max(abs(ip))
  scor[j]=ip[pred[j]]
  s[j]=sign(sum(xr[,pred[j]]*r))
 bhat[j]=ip[pred[j]]/sqrt(sum(xr[,pred[j]]^2))
}
return(list(pred=pred,s=s,scor=scor,bhat=bhat))
}


set.seed(1)
n <- 40
p <- 80
index <- sort(rep(1:(p/2), 2))
steps <- 10
sparsity <- 5
snr <- 3


    y <- rnorm(n)
    x <- matrix(rnorm(n*p), nrow=n)

  
      beta <- rep(0, p)
      beta[which(index %in% 1:sparsity)] <- snr
      y <- y + x %*% beta
  

    fit <- groupfs(x, y, index, maxsteps = steps)
    pvals <- groupfsInf(fit)


#test of size 1 groups
#

set.seed(1)
n <- 40
p <- 20
index <- sort(rep(1:(p/2), 2))
steps <- 10
sparsity <- 5
snr <- 3


 
    x <- matrix(rnorm(n*p), nrow=n)

  
      beta <- rep(0, p)
      beta[which(index %in% 1:sparsity)] <- snr
      y <- x %*% beta+sgma*rnorm(n)



fit <- groupfs(x, y, index=1:p, maxsteps = steps)

groupfsInf(fit)

xx=x
xx[,1:2]=0
xxx=scale_groups(xx,index)


fit2=myfs(x,y) #my old fs
fit3=fs(x,y,norm=FALSE)   #current
fit4=lars(x,y,type="step",norm=TRUE)

fit2$pred
fit3$act
fit4$act
max(abs(fit2$pred[1:38]-fit3$action[1:38]))
# They differ at the last entry, but that's OK (not well-defined when p>n)

max(abs(fit3$action[1:38]-unlist(fit4$action[1:38])))
# These don't always match, they make different selections at times. What
# is the lars function doing, in type="step"?

rbind(fit$act,fit2$pred[1:10],fit3$act[1:10])

fsInf(fit3,sigma=1)
fsInf(fit3,sigma=1,bits=200)

#minmodel=lm(y~1)
#step(minmodel,direction="forward")   #R step
#fm = step(minmodel, direction='forward', scope=(~x[,1]+x[,2]+x[,3]+x[,4]+x[,5]+x[,6]+x[,7]+x[,8]+x[,9]+x[,10]))
# fm$terms



      





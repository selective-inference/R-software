
library(intervals)
source("../selectiveInference/R/funs.common.R")
source("../selectiveInference/R/funs.groupfs.R")
source("../selectiveInference/R/funs.quadratic.R")
#source("../selectiveInference/R/funs.fs.R")
#source("../selectiveInference/R/funs.lar.R")
library(selectiveInference)#,lib.loc="/Users/tibs/dropbox/git/R/mylib")



myfs=function(x,y,nsteps=min(nrow(x),ncol(x))){
p=ncol(x)
# fs by minimizing scaled ip
# first center x and y
x=scale(x,T,F)
y=y-mean(y)
pred=s=scor=bhat=rep(NA,nsteps)
   ip=t(x)%*%y/sqrt(diag(t(x)%*%x))
  pred[1]=which.max(abs(ip))
 s[1]=sign(sum(x[,pred[1]]*y))
scor[1]=ip[pred[1]]
bhat[1]=ip[pred[1]]/sqrt(sum(x[,pred[1]]^2))

r=lsfit(x[,pred[1]],y)$res
for(j in 2:nsteps){
  mod=pred[1:(j-1)]
  r= lsfit(x[,mod],r)$res
  xr= lsfit(x[,mod],x)$res
  ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
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
  
 fit <- groupfs(x, y, index=1:p, maxsteps = steps)
   a=fs(x,y,normalize=T)

fit2=myfs(x,y) #my old fs
fit3=fs(x,y,norm=FALSE)   #current

rbind(fit$act,fit2$pred[1:10],fit3$act[1:10])
fsInf(fit3,sigma=1)
fsInf(fit3,sigma=1,bits=200)

#minmodel=lm(y~1)
#step(minmodel,direction="forward")   #R step
#fm = step(minmodel, direction='forward', scope=(~x[,1]+x[,2]+x[,3]+x[,4]+x[,5]+x[,6]+x[,7]+x[,8]+x[,9]+x[,10]))
# fm$terms

      





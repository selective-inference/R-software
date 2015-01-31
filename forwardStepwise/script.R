source("morefuns.R")


set.seed(33)
n=20
p=10
sigma=1

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

#generate y
beta=c(4,2,rep(0,p-2))
y=x%*%beta+sigma*rnorm(n)


a=myfs(x,y)

aa=myfs.pval(a,x,y,sigma)


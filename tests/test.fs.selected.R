library(selectiveInference)
library(lars)

set.seed(32)

n=50
p=10
sigma=1

x = as.matrix(read.table("X.csv", sep=',', header=FALSE))
Y = as.numeric(read.table("Y.csv", sep=',', header=FALSE)[,1])

beta=c(5,4,3,2,1,rep(0,p-5))
mu=x%*%beta

y=mu+Y
fsfit=fs(x,y,norm=TRUE, intercept=TRUE)
out = fsInf_maxZ(fsfit,sigma=sigma)



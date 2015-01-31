


library(lars)
library(covTest,lib.loc="lib.covTest")
library(larsInference,lib.loc="lib.larsInference")

# this uses an older version of genlasso (older than CRAN version)
library(genlasso,lib.loc="/Users/tibs/Dropbox/spacings/jasa-rev/mylib")

set.seed(33)
n=20
p=10
sigma=1

x=matrix(rnorm(n*p),n,p)
x=scale(x,T,T)/sqrt(n-1)

#generate y
beta=c(3,3,rep(0,p-2))
y=x%*%beta+sigma*rnorm(n)

aa=lars(x,y,norm=F,type="lasso")
aaa=covTest(aa,x,y,sigma=sigma)

aa2=lars(x,y,norm=F,type="lar")

b=larsInference(aa2,x,y,sigma=sigma)
pv.spacing[ii,1:nsteps]=b$pvalues.forw


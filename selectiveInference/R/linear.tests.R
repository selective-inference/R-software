# robs.test <- function() {
#   n <- 100
#   p <- 200
#   
#   set.seed(11332)
#   
#   y <- matrix(rnorm(n),ncol=1) # rand N(0,1) response
#   X <- matrix(rnorm(p*n),ncol = p) # p rand N(0,1) predictors
#   
#   X=scale(X,T,T)/sqrt(n-1)
#   lambda=1
#   sigma = estimateSigma(X,y)$sigmahat
#   
#   las <- glmnet(X,y,family="gaussian",alpha=1,standardize=F,intercept=T)
#   hbeta <- as.numeric(coef(las,x=X,y=y,s=lambda/n,exact=TRUE,intercept=T))
#   
#   
#   return(fixedLassoInf(X,y,hbeta[-1],lambda,family="gaussian",type="partial",intercept=T,sigma=sigma))
# }
# 
# 
# ## Tests partial inf for X and y randomly generated from N(0,1)
# nullTest <- function(X,y,lambda,intercept,type=c("full","partial")) {
#   n=nrow(X)
#   X=scale(X,T,T)/sqrt(n-1)
#   
#   sigma = estimateSigma(X,y)$sigmahat
#   
#   las <- glmnet(X,y,family="gaussian",alpha=1,standardize=F,intercept=intercept)
#   hbeta <- as.numeric(coef(las,x=X,y=y,s=lambda/n,exact=TRUE,intercept=intercept))
#   
#   if (type=="partial" || intercept==F) hbeta = hbeta[-1]
#   
#   return(fixedLassoInf(X,y,hbeta,lambda,family="gaussian",type=type,intercept=intercept,sigma=sigma))
# }
# 
# ## Test partial inf for X and y where 10 variables are y with random additive N(0,0.5) noise
# corrTest <- function(X,y,lambda,intercept,type=c("full","partial")) {
#   n=nrow(X)
#   corr.X = rep(y,10) + matrix(rnorm(n*10,0,0.5),ncol = 10)
#   X = cbind(corr.X,X)
#   X=scale(X,T,T)/sqrt(n-1)
#   
#   las <- glmnet(X,y,family="gaussian",alpha=1,standardize=F,intercept=intercept)
#   hbeta <- as.numeric(coef(las,x=X,y=y,s=lambda/n,exact=TRUE,intercept=intercept))
#   
#   sigma = estimateSigma(X,y)$sigmahat
#   
#   if (type=="partial" || intercept==F) hbeta = hbeta[-1]
#   
#   return(fixedLassoInf(X,y,hbeta,lambda,family="gaussian",type=type,intercept=intercept,sigma=sigma))
# }
# 
# ## QQ plot of p-values for all null data now that bug fix is implemented
# partial.qq.test <- function() {
#   n <- 100
#   p <- 200
#   
#   lambda=1
#   
#   null.int.pvs <- c()
#   corr.int.pvs <- c()
#   null.pvs <- c()
#   corr.pvs <- c()
#   for(i in 1:25) {
#     y <- matrix(rnorm(n),ncol=1) # rand N(0,1) response
#     X <- matrix(rnorm(p*n),ncol=p) # p rand N(0,1) predictors
#     
#     null <- nullTest(X,y,lambda,F,type="partial")
#     corr <- corrTest(X,y,lambda,F,type="partial")
#     null.pvs <- c(null.pvs,null$pv,recursive=T)
#     corr.pvs <- c(corr.pvs,corr$pv,recursive=T)
#     null.int <- nullTest(X,y,lambda,T,type="partial")
#     corr.int <- corrTest(X,y,lambda,T,type="partial")
#     null.int.pvs <- c(null.int.pvs,null.int$pv,recursive=T)
#     corr.int.pvs <- c(corr.int.pvs,corr.int$pv,recursive=T)
#   }
#   
#   qqplot(x=runif(length(null.pvs)),y=null.pvs,xlab="Expected",ylab="Observed",main="Partial Coef. Null X w/o Intercept")
#   abline(0,1)
#   qqplot(x=runif(length(corr.pvs)),y=corr.pvs,xlab="Expected",ylab="Observed",main="Partial Coef. 10 Corr. X w/o Intercept")
#   abline(0,1)
#   qqplot(x=runif(length(null.int.pvs)),y=null.int.pvs,xlab="Expected",ylab="Observed",main="Partial Coef. Null X w/ Intercept")
#   abline(0,1)
#   qqplot(x=runif(length(corr.int.pvs)),y=corr.int.pvs,xlab="Expected",ylab="Observed",main="Partial Coef. 10 Corr. X w/ Intercept")
#   abline(0,1)
# }
# 
# ## QQ plot of p-values for all null data now that bug fix is implemented
# pop.qq.test <- function() {
#   n <- 100
#   p <- 200
#   
#   lambda=1
#   
#   null.int.pvs <- c()
#   corr.int.pvs <- c()
#   null.pvs <- c()
#   corr.pvs <- c()
#   for(i in 1:25) {
#     y <- matrix(rnorm(n),ncol=1) # rand N(0,1) response
#     X <- matrix(rnorm(p*n),ncol=p) # p rand N(0,1) predictors
#     
#     null <- nullTest(X,y,lambda,F,type="full")
#     corr <- corrTest(X,y,lambda,F,type="full")
#     null.pvs <- c(null.pvs,null$pv,recursive=T)
#     corr.pvs <- c(corr.pvs,corr$pv,recursive=T)
#     null.int <- nullTest(X,y,lambda,T,type="full")
#     corr.int <- corrTest(X,y,lambda,T,type="full")
#     null.int.pvs <- c(null.int.pvs,null.int$pv,recursive=T)
#     corr.int.pvs <- c(corr.int.pvs,corr.int$pv,recursive=T)
#   }
#   
#   qqplot(x=runif(length(null.pvs)),y=null.pvs,xlab="Expected",ylab="Observed",main="Pop Coef. Null X w/o Intercept")
#   abline(0,1)
#   qqplot(x=runif(length(corr.pvs)),y=corr.pvs,xlab="Expected",ylab="Observed",main="Pop Coef. 10 Corr. X w/o Intercept")
#   abline(0,1)
#   qqplot(x=runif(length(null.int.pvs)),y=null.int.pvs,xlab="Expected",ylab="Observed",main="Pop Coef. Null X w/ Intercept")
#   abline(0,1)
#   qqplot(x=runif(length(corr.int.pvs)),y=corr.int.pvs,xlab="Expected",ylab="Observed",main="Pop Coef. 10 Corr. X w/ Intercept")
#   abline(0,1)
# }
# 
# 
# 
# 
# ## QQ plot of p-values for data with correlated x now that bug fix implemented
# power.partial.pval.dist <- function(n,p,intercept=T,lambda=1) {
#   pvs <- c()
#   for(i in 1:10) {
#     a <- powerPartialTest(n,p,intercept,lambda)
#     ps <- a$pv
#     pvs <- c(pvs,ps,recursive=T)
#   }
#   qqplot(x=runif(length(pvs)),y=pvs,xlab="Expected",ylab="Observed",main="Partial Coef. 10 Corr. X")
#   abline(0,1)
# }
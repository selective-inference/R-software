\name{selectiveInference}
\alias{selectiveInference}
\docType{package}
\title{
Tools for selective inference
}
\description{
Functions to perform post-selection inference for forward
stepwise regression, least angle regression, the lasso and the 
many normal means problem
}
\details{
\tabular{ll}{
Package: \tab selectiveInference\cr
Type: \tab Package\cr
License: \tab GPL-2\cr
}

This package provides tools for inference after selection, in forward stepwise 
regression, least angle regression, the lasso, and the many normal means problem. 
The functions compute p-values and selection intervals that properly account for 
the inherent selection carried out by the procedure. These have exact finite sample 
type I error and coverage under Gaussian errors.

This R package was developed as part of the selective inference software project
in Python and R:

\url{https://github.com/selective-inference}

Some of the R code in this work is a modification of Python code from this
repository. Here is the current selective inference software team:

Yuval Benjamini, 
Leonard Blier, 
Will Fithian, 
Jason Lee, 
Joshua Loftus,
Joshua Loftus, Stephen Reid,
Dennis Sun,
Yuekai Sun,
Jonathan Taylor,
Xiaoying Tian,
Ryan Tibshirani,
Rob Tibshirani

The main functions included in the package are:
\code{\link{fs}},
\code{\link{fsInf}},
\code{\link{lar}},
\code{\link{larInf}},
\code{\link{fixedLassoInf}},
\code{\link{manyMeans}}
}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid

Maintainer: Rob Tibshirani <tibs@stanford.edu>
}

\references{
Ryan Tibshirani, Jonathan Taylor, Richard Lockhart, Rob Tibshirani (2014).
Exact post-selection inference for sequential regression procedures. arXiv:1401.3889.

Jason Lee, Dennis Sun, Yuekai Sun, Jonathan Taylor (2013). 
Exact post-selection inference, with application to the lasso. arXiv:1311.6238.

Stephen Reid, Jonathan Taylor, Rob Tibshirani (2014).
Post-selection point and interval estimation of signal sizes in Gaussian samples.
arXiv:1405.3340.
}

\examples{
#NOT RUN
# forward stepwise:
#
#generate some data
#set.seed(43)
#n=50
#p=10
#sigma=.7
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,T)/sqrt(n-1)
#beta=c(4,2,0,0,rep(0,p-4))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#
#first run forward stepwise
# fsfit=fs(x,y)
#
# forward stepwise inference for each successive entry of a predictor; 
#   sigma estimated from mean squared residual
#   
# aa=fsInf(fsfit)
##
#  lasso with fixed lambda
#
#set.seed(43)
#n=50
#p=10
#sigma=.7
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,T)/sqrt(n-1)
#beta=c(4,2,0,0,rep(0,p-4))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#
# first run  glmnet
#gfit=glmnet(x,y,standardize=F)
#lam = .1
#extract coef for a given lam; Note the 1/n factor in s!
#bhat = coef(gfit, s=lam/n, exact=TRUE)[-1]

# compute fixed lambda p-values and selection intervals
#aa=fixedLassoInf(x,y,bhat,lam,sigma=sigma)
#
##least angle regression from mean squared residual
#set.seed(43)
#n=50
#p=10
#sigma=.7
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,T)/sqrt(n-1)
#beta=c(3,2,0,0,rep(0,p-4))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#first run lar
# larfit=lar(x,y)
#
#lar inference for each successive entry of a predictor; sigma estimated
#  from mean squared residual from least squares fit
# aa=larInf(larfit)

##
##many normal means

#set.seed(12345)
#n = 100 # sample size
#signal = 3 # signal size
#mu = c(rep(signal, floor (n/5)), rep(0, n-floor(n/5))) # 20% of elements get the signal; rest 0
#y = mu + rnorm (n, 0, 1)
#mmObj = manyMeans(y, bh.q=0.1)
#
#
# Forward stepwise with grouped variables
#set.seed(1)
#n <- 40
#p <- 20
#index <- sort(rep(1:(p/2), 2))
#steps <- 10
#sparsity <- 5
#snr <- 3
#   x <- matrix(rnorm(n*p), nrow=n)
#      beta <- rep(0, p)
#      beta[which(index %in% 1:sparsity)] <- snr
#      y <-  x %*% beta+rnorm(n)
#fit <- groupfs(x, y, index=1:p, maxsteps = steps)
#out<- groupfsInf(fit)



##estimation of sigma for use in fsInf or larInf or fixedLassoInf
#
#set.seed(43)
#n=50
#p=10
#sigma=.7
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,F)
#beta=c(3,2,0,0,rep(0,p-4))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)

#out=estimateSigma(x,y)

##estimation of lambda  for use in fixedLassoInf
set.seed(43)
#n=50
#p=10
#sigma=.7
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,F)
#beta=c(3,2,0,0,rep(0,p-4))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#
#estimate lambda using the  known value of sigma
#lamhat=estimateLambda(x,sigma=.7)
#first estimate  sigma
#sigmahat=estimateSigma(x,y)$sigmahat
#lamhat=estimateLambda(x,sigma=sigmahat)
#

}
\keyword{ package }



\name{selectiveInference}
\alias{selectiveInference}
\docType{package}
\title{
Tools for selective inference
}
\description{
Functions to perform inference after selection for forward
stepwise regression, least angle regression, the lasso and the 
many normal means problem
}
\details{
\tabular{ll}{
Package: \tab selectiveInference\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-09-01\cr
License: \tab GPL-2\cr
}

This package provides tools for inference after selection, in forward stepwise 
regression, least angle regression, the lasso, and the many normal means problem. 
The functions compute p-values and selection intervals that properly account for 
the inherent selection carried out by the procedure. These have exact finite sample 
type I error and coverage under Gaussian errors.

This R code was developed in coordination with an analogous Python package developed 
by Jonathan Taylor and others, see http://github.com/selective-inference/selective-inference.

The main functions included in the package are:
\code{\link{fs}},
\code{\link{fsInf}},
\code{\link{lar}},
\code{\link{larInf}},
\code{\link{fixedLassoInf}},
\code{\link{manyMeans}}
}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Stephen Reid

Maintainer: Rob Tibshirani <tibs@stanford.edu>
}

\references{
Jonathan Taylor, Richard Lockhart, Ryan Tibshirani, Rob Tibshirani (2014).
Exact post-selection inference for sequential regression procedures. arXiv:1401.3889.

Jason Lee, Dennis Sun, Yuekai Sun, Jonathan Taylor (2014). 
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

}
\keyword{ package }



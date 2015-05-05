\name{selectiveInference}
\alias{selectiveInference}
\docType{package}
\title{
Functions for Selective Inference
}
\description{
A collection of functions that performs inference after selection for forward stepwise regression,
the lasso, least angle regression and the many normal means problem.
}
\details{
\tabular{ll}{
Package: \tab selectiveInference\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-01-05\cr
License: \tab GPL-2\cr
}
Ths collection of functions  performs inference after selection for forward stepwise regression,
the lasso, least angle regression and the many normal means problem.
The functions compute p-values and selection intervals that properly account for the inherence slection
carried out by the procedure. These have exact finite sample  type I error and coverage
under Gaussian errors.

The following functions are included in the package:

\code{\link{forwardStep}},
\code{\link{forwardStepInf}},
\code{\link{fixedLassoInf}},
\code{\link{lar}},
\code{\link{predict.lar}},
\code{\link{larInf}},
\code{\link{forwardStop}},
\code{\link{manyMeans}}

}

\author{Joshua Loftus, Stephen Reid, Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani


Maintainer: Rob Tibshirani <tibs@stanford.edu>
}

\references{
Jason D. Lee, Dennis L. Sun, Yuekai Sun, Jonathan E. Taylor (2014). Exact post-selection inference, with application to the lasso. arXiv:1311.6238

Stephen Reid, Jonathan Taylor, Robert Tibshirani (2014).
Post-selection point and interval estimation of signal sizes in Gaussian samples.
arXiv:1405.3340

Jonathan Taylor, Richard Lockhart, Ryan Tibshirani, Rob Tibshirani (2014).
Exact Post-selection Inference for Forward Stepwise and Least Angle Regression. arXiv:1401.3889
}
\examples{
#NOT RUN
# forward stepwise:
#set.seed(43)
#n=50
#p=10
#sigma=.7
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,F)
#beta=c(3,2,0,0,rep(0,p-4))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#first run forward stepwise
# fsfit=forwardStep(x,y)
#
# forward stepwise inference for each successive entry of a predictor; sigma estimated
#  from mean squared residual
# aa=forwardStepInf(fsfit,x,y)
##
#  lasso with fixed lambda
#
#set.seed(43)
#n=50
#p=10
#sigma=.7
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,T)/sqrt(n-1)
#beta=c(3,2,0,0,rep(0,p-4))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
# first run  glmnet
#gfit=glmnet(x,y,standardize=F)
#lam = .1
#extract coef for a given lam; Note the 1/n factor!
#bhat = coef(gfit, s=lam/n, exact=TRUE)[-1]

# compute fixed lambda p-values and selection intervals
#aa=fixedLassoInf(x,y,bhat,lam,sigma)

##least angle regression
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
# aa=larInf(larfit,x,y)

##
##many normal means

#set.seed(12345)
#n = 100 # sample size
#signal = 3 # signal size
#mu = c(rep(signal, floor (n/5)), rep(0, n-floor(n/5))) # 20% of elements get the signal; rest 0
#y = mu + rnorm (n, 0, 1)
#mmObj = manyMeans(y, bh.q=0.1)


##estimation of sigma
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



\name{estimateLambda}
\alias{estimateLambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates the lasso tuning parameter lambda.
}
\description{
Estimates the lasso tuning parameter lambda, for use in the selectiveInference 
package
}
\usage{
estimateLambda(x, sigma, nsamp=1000)
}
\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{sigma}{
Estimate of error standard deviation
}
\item{nsamp}{Number of Monte carlo samples used for the estimation.}
}
\details{
This function estimates the lasso tuning parameter lambda, using the estimate
2*E(||X^T eps||_infty)  where eps ~ N(0,sigma^2), a vector of length n.
This estimate was proposed by Negahban et  al (2012).
}
\value{
\item{sigmahat}{The estimate of sigma}
\item{df}{The degrees of freedom of lasso fit used}
}
\references{
Negahban, S. N.,
Ravikumar, P.,
Wainwright, M. J.
and Yu, B.
(2012).  A unified
framework for high-dimensional analysis of
M-estimators with decomposable regularizers.
Statistical Science vol. 27, p 538-557.
}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid}

\examples{
#NOT RUN
#set.seed(43)
#n=50
#p=10
#sigma=.7
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,F)
#beta=c(3,2,0,0,rep(0,p-4))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#
#estimate lambda usingthe  known value of sigma
#lamhat=estimateLambda(x,sigma=.7)
#
#first estimate  sigma
#sigmahat=estimateSigma(x,y)$sigmahat
#lamhat=estimateLambda(x,sigma=sigmahat)

#compare to estimate from cv

#out=cv.glmnet(x,y)
#out$lambda.min*n    #remember that value from glmnet must be
                     # multiplied by n, to make it comparable.
}


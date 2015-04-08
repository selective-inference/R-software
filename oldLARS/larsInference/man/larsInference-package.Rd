\name{larsInference-package}
\alias{larsInference-package}
\docType{package}
\title{
Compute p-values and confidence intervals after each step of the LAR
algorithm
}
\description{
 This package computes p-values and confidence intervals after each
     steps of the least angle regresson (LAR) algorithm.
}
\details{
\tabular{ll}{
Package: \tab larsInference\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-01-08\cr
License: \tab GPL-2\cr
}
Very simple to use. Takes output from a call to  \code{lars}.
Requires \code{lars} and \code{glasso} packages.
are included in this package.
Functions are:
\code{larsInference}
}
\author{Rob Tibshirani tibs@stanford.edu
}
\references{
A significance test for the lasso (2014). Lockhart, R.,  Taylor, J.,  Tibshirani (Ryan) and  Tibshirani (Robert)

Adaptive inference for Least Angle Regression and the  Lasso (2014). Taylor, J.  Lockhart, R.,   Tibshirani (Ryan) and  Tibshirani (Robert). 
}
\keyword{ LARS }
\seealso{larsInference}
\examples{

x=matrix(rnorm(100*10),ncol=10)
x=scale(x,TRUE,TRUE)/sqrt(99)

beta=c(4,rep(0,9))
y=x\%*\%beta+.4*rnorm(100)

a=lars(x,y)
larsInference(a,x,y)

}


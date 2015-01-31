\name{larsInference}
\alias{larsInference}
\title{Compute p-values and confidence intervals after each step of  the Least Angle Regression LAR algorithm.
}
\description{
This package computes p-values and confidence intervals after each step of  the least angle regresson (LAR) algorithm.
}
\usage{
larsInference(larsfit, x, y, sigma=NULL, nsteps=ncol(x),alpha=0.1,trace=FALSE)
}
\arguments{
  \item{larsfit}{
Result of a call to \code{lars} 
}
  \item{x}{
N by p matrix of predictors
     }
  \item{y}{
N-vector of outcome values
}
  \item{sigma}{
Estimate of the standard deviation of the noise. If a numerical value, that value if used.
Otherwise it is estimated from the MSE of CV-estimated model.
}
\item{nsteps}{
Total number of LAR steps to be analyzed.
}
\item{alpha}{
The desired (two-sided) confidence level
}
\item{trace}{
flag to provide info along the way
}

}
\details{
This package computes p-values and confidence intervals after each steps of  the least angle regresson (LAR) algorithm. It returns both exact pvalues and selection intervals based on the ``spacings'' method of Taylor et al (2014), and
also asymptotic p-values based on the covariance test of Lockhart et al (2014).
}
\value{
\item{pvalues.forward}{Exact (finite sample) pvalues for each predictor as it enters the model} 
\item{pvalues.remove}{Exact (finite sample) pvalues for each predictor when it is removed from the current the model}
\item{pvalues.forward.cov.unc}{Asymptotic  pvalues for each predictor as it enters the model,
 based on the covariance test of Lockhart et al (2014), using Exp(1) null at each step} 
\item{pvalues.forward.cov.cor}{Asymptotic  pvalues for each predictor as it enters the model,
based on the corrected covariance test of Taylor et al (2014)}
\item{selection.intervals}{ alpha-level confidence intervals for  each predictor at each LAR step.
   See Taylor et al (2014) for precise interpretation of these intervals}
\item{alpha}{The desired (two-sided) confidence level used}
  \item{sigma}{Estimate of sigma used}
}
\references{
A significance test for the lasso (2013). Lockhart, R.,  Taylor, J.,  Tibshirani (Ryan) and  Tibshirani (Robert)

Adaptive inference for Least Angle Regression and the  Lasso (2014). Taylor, J.  Lockhart, R.,   Tibshirani (Ryan) and  Tibshirani (Robert).

}
\author{
Rob Tibshirani
}
\note{
This function requires  \code{lars} and \code{genlasso} functions from the CRAN R collection.
}


\seealso{
lars
}
\examples{
set.seed(1234)
#generate data
x=matrix(rnorm(100*10),ncol=10)
#standardize
mx=colMeans(x); sx=sqrt(apply(x,2,var)); x=scale(x,mx,sx)/sqrt(9)
#generate outcome
beta=c(3,-2,rep(0,8))
y=x\%*\%beta+2*rnorm(100)
#run lars
a=lars(x,y,type="lar")

# do inference
larsInference(a,x,y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{least angle regression} 
\keyword{lasso}
\keyword{adaptive fitting} 

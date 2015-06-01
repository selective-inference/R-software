\name{larInf}
\alias{larInf}
\title{
Selective inference for least angle regression (LAR)
}
\description{
Compute p-values and selection intervals for least angle regression
}
\usage{
larInf(x, y, larfit, sigma=NULL, compute.si = TRUE, alpha = 0.1, one.sided = TRUE, nsteps = min(nrow(x), ncol(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
  \item{larfit}{
Object returned by \code{lar} function (not the \code{lars} function)!
}
  \item{sigma}{
Estimate of error standard deviation
}
  \item{compute.si}{
Compute selection intervals? Default TRUE. Can slow down the computation, if only p-values are required.
}
  \item{alpha}{
Significance level for  selection intervals (target is miscoverage alpha/2 in each tail)
}
  \item{one.sided}{
One-sided p-values? Default TRUE.
}
\item{nsteps}{
Maximum number of steps to take (default=min(n,p))
}
}
\details{
This function computes selective p-values and confidence intervals (selection intervals) for least angle regression.   The selection interval estimation involves numerical search and can be fragile:
if the observed statistics are too close to either end of the truncation interavls (Vm, Vp) (see ref), the interval of desired coverage
cannot be computed and is set to (-Inf, +Inf). The output \code{tailarea} gives the achieved Gaussian tail area
for reported  interval--- these should be close to alpha/2.}

\value{
  \item{pv}{p-values for each predictor}
\item{ci}{Selection intervals}
\item{tailarea}{Realized tail area for each  selection interval}
  \item{vm}{Lower bounds for functionals}
 \item{vp}{Upper bounds for functionals}
\item{sigma}{Value of error standard deviation (sigma) used}
\item{alpha}{Desired coverage (alpha/2 in each tail)}
\item{act}{Predictors entered at each step of LAR}
\item{pv.spacing}{p-values from asymptotic spacing test}
\item{pv.cov}{p-values from covariance (exponential) rtest}
\item{ForwardStopHat}{Estimated stopping point from forwardStop rule, approximately controlling FDR at alpha}
\item{call}{The call to forwardStepInf}
}

\references{
Jonathan Taylor, Richard Lockhart, Ryan Tibshirani, Rob Tibshirani (2014).
Exact Post-selection Inference for Forward Stepwise and Least Angle Regression. arXiv:1401.3889

}
\author{
Joshua Loftus, Stephen Reid, Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{lar}}
}
\examples{

#NOT RUN
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
# aa=larInf(x,y,larfit)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

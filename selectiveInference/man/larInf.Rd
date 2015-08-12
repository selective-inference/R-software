\name{larInf}
\alias{larInf}
\title{
Selective inference for least angle regression (LAR)
}
\description{
Compute p-values and selection intervals for least angle regression
}
\usage{
larInf(obj, sigma=NULL, alpha=0.1, k=NULL, type=c("active","all","aic"), 
                   gridfac=25, gridpts=1000, mult=2, ntimes=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{obj}{
Object returned by \code{lar} function (not the \code{lars} function)!
}
  \item{sigma}{Estimate of error standard deviation}
  \item{alpha}{
Significance level for  selection intervals (target is miscoverage alpha/2 in each tail). Eg alpha=.10 or .05.
}
  \item{k}{
See "type" argument below. Default NULL
}
\item{type}{Type of analysis desired: with "active" (default), p-values and confidence intervals are computed for each predictor as it is entered; with "All",  p-values and confidence intervals are computed
for model constructed after "k" steps; with "aic", the number of steps k is estimated the AIC
criterion, and then p-values and confidence intervals are computed
for model constructed after this number of steps. Note that BIC or other criteria can be specified 
via the argument "mult" below.
}
\item{gridfac}{Factor that multiplies the standard deviation, in determining the width of the grid we use for confidence interval construction. Default 25.}
\item{gridpts}{Number of grid points in the grid used  for confidence interval construction. Default 1000.}
\item{mult}{Multiplier for AIC-style penallty; default of 2 gives AIC: mult=log(nrow(x)) gives BIC}
\item{ntimes}{Number of steps for which the AIC-style criterion has to increase before minimizing point is declared. Default 2.}
}
\details{
This function computes selective p-values and confidence intervals (selection intervals) for least angle regression.   The selection interval estimation involves numerical search and can be fragile:
if the observed statistics are too close to either end of the truncation interavls (Vm, Vp) (see ref), the interval of desired coverage
cannot be computed and is set to (-Inf, +Inf). The output \code{tailarea} gives the achieved Gaussian tail area
for reported  interval--- these should be close to alpha/2.}

\value{

\item{type}{Type of analysis specified in call}
  \item{k}{Value of k specified in call}
\item{khat}{Estimated number of steps, if type="aic" was specified}
\item{pv}{Estimated pvalues}
 \item{ci}{Selection intervals}
\item{tailarea}{Realized tail area (lower and upper) for each  selection interval}
 \item{vlo}{Lower bounds for functionals}
 \item{vup}{Upper bounds for functionals}
\item{vmat}{Matrix of linear functionals used for inference}
\item{y}{Vector of outcomes}
\item{pv.spacing}{p-values from the spacing test (where we use M^+ as V^{lo})}
 \item{pv.asymp}{p-values from the asymptotic form of the spacing test (where we use \lambda_{k-1} (previous knot) in place of M^+)}
\item{pv.cov}{p-values from covariance test. The reference distribution is Exp(1) if sigma is unknown; otherwise it is F_(2,n-p)}
\item{vars}{Predictors entered at each step of forward stepwise}
\item{sign}{Signs of coefficients for predictors entered at each step of forward stepwise}
\item{sigma}{Value of error standard deviation (sigma) used}
\item{alpha}{Desired coverage (alpha/2 in each tail)}
\item{call}{The call to larInf}
}

\references{
Jonathan Taylor, Richard Lockhart, Ryan Tibshirani, Rob Tibshirani (2014).
Exact Post-selection Inference for Forward Stepwise and Least Angle Regression. arXiv:1401.3889

}
\author{Ryan Tibshirani, Rob Tibshirani,  Jonathan Taylor, Max G'Sell, Stephen Reid}



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
# aa=larInf(larfit)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

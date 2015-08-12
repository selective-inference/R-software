\name{fsInf}
\alias{fsInf}
\title{
Selective inference for forward stepwise regression
}
\description{
Computes exact p-values and selection intervals for forward stepwise regression
}
\usage{
fsInf(obj, sigma=NULL,  alpha = 0.1, k=NULL, type=c("active","all","aic"), 
                  gridfac=25, gridpts=1000, mult=2, ntimes=2) 
}
\arguments{
  \item{obj}{
Object returned by fs
}

\item{sigma}{
Estimate of error standard deviation; default is value returned by fs.
}
 
  \item{alpha}{
Significance level for  selection intervals (target is miscoverage alpha/2 in each tail). Eg alpha=.10 or .05.
}
  \item{k}{
See "type" argument below. Default NULL;
}
\item{type}{Type of analysis desired: with "active" (default), p-values and confidence intervals are computed for each predictor as it is entered; with "All",  p-values and confidence intervals are computed
for model constructed after "k" steps; with "aic", the number of steps k is estimated the AIC
criterion, and then p-values and confidence intervals are computed
for model constructed after this number of steps. Note that BIC or other criteria can be specified 
htought the argument "mult" below.
}
\item{gridfac}{Factor that multiplies the standard deviation, in determining the width of the grid we use for confidence interval construction. Default 25.}
\item{gridpts}{Number of grid points in the grid used  for confidence interval construction. Default 1000.}
\item{mult}{Multiplier for AIC-style penallty; default of 2 gives AIC: mult=log(nrow(x)) gives BIC}
\item{ntimes}{Number of steps for which AIC-style criterion has to increase before minimizing point is declared}
}
\details{
This function computes selective p-values and confidence intervals (selection intervals) for forward stepwise regression. The default is to report the results for each predictor after its entry into the model.
See the "type" argument for other options.
  The selection interval estimation involves numerical search and can be fragile:
if the observed statistics are too close to either end of the truncation interavls (Vm, Vp) (see ref), the interval of desired coverage
cannot be computed and is set to (-Inf, +Inf). The output \code{tailarea} gives the achieved Gaussian tail area
for reported  interval--- these should be close to alpha/2.
}

\value{
\item{type}{Type of analysis specified in call}
  \item{k}{Value of k specified in call}
\item{khat}{Estimated number of steps, if type="aic" was used}
\item{pv}{Estimated pvalues}
 \item{ci}{Selection intervals}
\item{tailarea}{Realized tail area for each  selection interval}
 \item{vlo}{Lower bounds for functionals}
 \item{vup}{Upper bounds for functionals}
\item{vmat}{Matrix of linear functionals used for inference}
\item{y}{Vector of outcomes}
\item{vars}{Predictors entered at each step of forward stepwise}
\item{sign}{Signs of coefficients for predictors entered at each step of forward stepwise}
\item{which.steps }{Which FS steps were analyzed}
\item{step.ind}{Step indicator for each row of constraint matrix A (returned below)}
\item{ForwardStopHat}{Estimated stopping point from forwardStop rule, approximately controlling FDR at alpha}
\item{sigma}{Value of error standard devarion (sigma) used}
\item{alpha}{Desired coverage (alpha/2 in each tail)}
\item{call}{The call to fsInf}
}
\references{
Jonathan Taylor, Richard Lockhart, Ryan Tibshirani, Rob Tibshirani (2014).
Exact Post-selection Inference for Forward Stepwise and Least Angle Regression. arXiv:1401.3889


Joshua R. Loftus, Jonathan E. Taylor (2014).
A significance test for forward stepwise
model selection. arXiv:1405.3920v1
}

\author{Ryan Tibshirani, Rob Tibshirani,  Jonathan Taylor, Max G'Sell, Stephen Reid}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{fs}}
}
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
#first run forward stepwise
# fsfit=fs(x,y)
#
# forward stepwise inference for each successive entry of a predictor; sigma estimated
#  from mean squared residual
# aa=fsInf(fsfit)

# forward stepwise inference for fixed model of size 4, testing deletion of each predictor; 
# known value of sigma used
# aa2=fsInf(fsfit,x,y,sigma=sigma,type="all",k=4)

# forward stepwise inference for  model chosen by AIC, testing deletion of each predictor; 
# known value of sigma used
# aa3=fsInf(fsfit,sigma=sigma,type="aic")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

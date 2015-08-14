\name{fsInf}
\alias{fsInf}
\title{
Selective inference for forward stepwise regression
}
\description{
Computes p-values and confidence intervals for forward 
stepwise regression
}
\usage{
fsInf(obj, sigma=NULL, alpha=0.1, k=NULL, type=c("active","all","aic"), 
      gridrange=c(-100,100), gridpts=1000, mult=2, ntimes=2, verbose=FALSE) 
}
\arguments{
  \item{obj}{
Object returned by \code{\link{fs}} function
}
\item{sigma}{
Estimate of error standard deviation. If NULL (default), this is estimated 
using the mean squared residual of the full least squares fit when n >= 2p, and 
using the standard deviation of y when n < 2p. In the latter case, the user 
should use \code{\link{estimateSigma}} function for a more accurate estimate
}
\item{alpha}{
Significance level for confidence intervals (target is miscoverage alpha/2 in each tail)
}
\item{k}{
See "type" argument below. Default is NULL, in which case k is taken to be the
the number of steps computed in the forward stepwise path
}
\item{type}{
Type of analysis desired: with "active" (default), p-values and confidence intervals are 
computed for each predictor as it is entered into the active step, all the way through
k steps; with "all", p-values and confidence intervals are computed for all variables in 
the active model after k steps; with "aic", the number of steps k is first estimated using 
a modified AIC criterion, and then the same type of analysis as in "all" is carried out for 
this particular value of k.  

Note that the AIC scheme is defined to choose a number of steps k after which the AIC criterion 
increases \code{ntimes} in a row, where \code{ntimes} can be specified by the user (see below). 
Under this definition, the AIC selection event is characterizable as a polyhedral set, and hence
the extra conditioning can be taken into account exactly. Also note that an analogous BIC scheme
can be specified through the \code{mult} argument (see below)
}
\item{gridrange}{
Grid range for constructing confidence intervals, on the standardized scale
}
\item{gridpts}{
Number of grid points for constructing confidence intervals, within the specified
grid range}
\item{mult}{Multiplier for the AIC-style penalty. Hence a value of 2 (default) 
gives AIC, whereas a value of log(n) would give BIC}
\item{ntimes}{Number of steps for which AIC-style criterion has to increase before 
minimizing point is declared}
\item{verbose}{Print out progress along the way? Default is FALSE}
}

\details{
This function computes selective p-values and confidence intervals (selection intervals)
for forward stepwise regression. The default is to report the results for 
each predictor after its entry into the model. See the "type" argument for other options.
The confidence interval construction involves numerical search and can be fragile:
if the observed statistics are too close to either end of the truncation intervals 
(vlo and vup, see references), then the interval of desired coverage
cannot be computed and is set to [-Inf, +Inf]. The output \code{tailarea} gives 
the achieved Gaussian tail area for reported interval---these should be close to alpha/2,
and can be used for error-checking purposes.
}

\value{  
\item{type}{Type of analysis (active, all, or aic)}
\item{k}{Value of k specified in call}
\item{khat}{When type is "active", this is an estimated stopping point
declared by \code{\link{forwardStop}}; when type is "aic", this is the
value chosen by the modified AIC scheme}
\item{pv}{P-values for active variables}
\item{ci}{Confidence intervals}
\item{tailarea}{Realized tail areas (lower and upper) for each confidence interval}
\item{vlo}{Lower truncation limits for statistics}
\item{vup}{Upper truncation limits for statistics}
\item{vmat}{Linear contrasts that define the observed statistics}
\item{y}{Vector of outcomes}
\item{vars}{Variables in active set}
\item{sign}{Signs of active coefficients}
\item{alpha}{Desired coverage (alpha/2 in each tail)}
\item{sigma}{Value of error standard deviation (sigma) used}
\item{call}{The call to fsInf}
}

\references{
Jonathan Taylor, Richard Lockhart, Ryan Tibshirani, Rob Tibshirani (2014).
Exact post-selection inference for sequential regression procedures. arXiv:1401.3889.

Joshua Loftus and Jonathan Taylor (2014). A significance test for forward stepwise
model selection. arXiv:1405.3920.
}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Stephen Reid}

\seealso{\code{\link{fs}}}

\examples{
#NOT RUN
#set.seed(433)

#n=50
#p=10
#sigma=.7
#x=matrix(rnorm(n*p),n,p)

#x=scale(x,T,T)/sqrt(n-1)
#beta=c(4,2,0,0,rep(0,p-4))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#first run forward stepwise
# fsfit=fs(x,y)
#
# forward stepwise inference for each successive entry of a predictor; 
#
#  sigma estimated from mean squared residual
#  
# aa=fsInf(fsfit)

# forward stepwise inference for fixed model of size 4, testing deletion of each predictor; 
# known value of sigma used
# aa2=fsInf(fsfit,sigma=sigma,type="all",k=4)

# forward stepwise inference for  model chosen by AIC, testing deletion of each predictor; 
# known value of sigma used
# aa3=fsInf(fsfit,sigma=sigma,type="aic")
#
# states data example: predict life expectancy
#x=state.x77[,-4]
#y=state.x77[,4]
#x=scale(x,T,T)
#n=nrow(x)
#  estimate sigma from cross-validated lasso fit
#cvf=cv.glmnet(x,y)
#sigmahat=estimateSigma(x,y,stand=F)$sigmahat
#
# run forward stepwise
#fsfit=fs(x,y,normalize=F)
#
# compute p-values and confidence intervals
# fsInf(fsfit,sigma=sigmahat)
}
\name{forwardStepInf}
\alias{forwardStepInf}
\title{
Selective inference for forward stepwise regression
}
\description{
Computes exact p-values and selection intervals for forward stepwise regression
}
\usage{
forwardStepInf(fsfit, x, y, sigma=NULL, nsteps = NULL, alpha = 0.1, fixed.step = NULL, aic.stop = FALSE, trace = F, compute.si = TRUE, one.sided = TRUE)
}
\arguments{
  \item{fsfit}{
Object returned by forwardStep
}
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
\item{sigma}{
Estimate of error standard deviation; default is value returned by forwardStep.
}
  \item{nsteps}{
Number of steps to analyze. Default is min(k,20) where k is the  number used in forwardStep.

}
  \item{alpha}{
Significance level for  selection intervals (target is miscoverage alpha/2 in each tail)
}
  \item{fixed.step}{
If an integer k, results for model after step k are computed. If NULL and aic.stop is FALSE (default),
results for each FS step are computed. 
}
  \item{aic.stop}{
 Use fixed step prescibed by  AIC? Default FALSE
}
  \item{trace}{
Print out more details along the way.
}
  \item{compute.si}{
Compute selection intervals? Default TRUE. Can slow down the computation, if only p-values are required.
}
  \item{one.sided}{
One-sided p-values? Default TRUE.
}
}
\details{
This function computes selective p-values and confidence intervals (selection intervals) for forward stepwise regression. The default is to report the results for each predictor after its entry into the model.
If fixed.step is not null  or aic.stop=T, the results are for the removal  of each predictor
from the fixed model specified.  The selection interval estimation involves numerical search and can be fragile:
if the observed statistics are too close to either end of the truncation interavls (Vm, Vp) (see ref), the interval of desired coverage
cannot be computed and is set to (-Inf, +Inf). The output \code{tailarea} gives the achieved Gaussian tail area
for reported  interval--- these should be close to alpha/2.
}

\value{
  \item{pv}{p-values for each predictor}
  \item{vm}{Lower bounds for functionals}
 \item{vp}{Upper bounds for functionals}
\item{ci}{Selection intervals}
\item{tailarea}{Realized tail area for each  selection interval}
\item{pred}{Predictors entered at each step of forward stepwise}
\item{which.steps }{Which FS steps were analyzed}
\item{step.ind}{Step indicator for each row of constraint matrix A (returned below)}
\item{ForwardStopHat}{Estimated stopping point from forwardStop rule, approximately controlling FDR at alpha}
\item{alpha}{Desired coverage (alpha/2 in each tail)}
\item{sigma}{Value of error standard devarion (sigma) used}
\item{one.sided}{Are the pvalues one-sided}
\item{A}{Constraint matrix}
\item{b}{Right hand side vector for constraints}
\item{call}{The call to forwardStepInf}
}
\references{
Jonathan Taylor, Richard Lockhart, Ryan Tibshirani, Rob Tibshirani (2014).
Exact Post-selection Inference for Forward Stepwise and Least Angle Regression. arXiv:1401.3889


Joshua R. Loftus, Jonathan E. Taylor (2014).
A significance test for forward stepwise
model selection. arXiv:1405.3920v1
}
\author{
Joshua Loftus, Stephen Reid,  Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{forwardStep}}
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
# fsfit=forwardStep(x,y)
#
# forward stepwise inference for each successive entry of a predictor; sigma estimated
#  from mean squared residual
# aa=forwardStepInf(fsfit,x,y)

# forward stepwise inference for fixed model of size 4, testing deletion of each predictor; 
# known value of sigma used
# aa2=forwardStepInf(fsfit,x,y,sigma=sigma,fixed.step=4)

# forward stepwise inference for  model chosen by AIC, testing deletion of each predictor; 
# known value of sigma used
# aa3=forwardStepInf(fsfit,x,y,sigma=sigma, aic.stop=T)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

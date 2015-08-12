\name{fixedLassoInf}
\alias{fixedLassoInf}

\title{
Inference for lasso with fixed lambda 
}
\description{
Compute p-values and selection intervals for lasso fit with fixed value of regularization parameter
lambda
}
\usage{
fixedLassoInf(x, y, bhat, lambda, coeftype = c("partial", "full"), sigma=NULL, alpha = 0.1, verbose = FALSE, compute.si = TRUE, tol.beta = 1e-5,
tol.kkt=0.05, one.sided = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
  \item{bhat}{
Estimated lasso coefs (eg from glmnet). Intercept is not included.

 Be careful!
 This function uses the "standard" lasso objective
(1/2)*RSS +lambda*sum(abs(beta)).  In contrast, glmnet uses a factor (1/2n) in the first term.
So after running glmnet, to extract the bhat corresponding to a value lambda, you need to use bhat=coef(fit, s=lambda/n)[-1],
where fit is the object returned by glmnet ([-1] removes the intercept). 
}
  \item{lambda}{
Value of lambda used to compute bhat
}
\item{coeftype}{Parameter type for selection intervals: default is "partial"--- the partial regression coefficients in the active set of predictors;
the other option is "full"-- the population regression coefficients, assuming the model with all predictors in correct. The latter is not allowed when p>n.}
  \item{sigma}{
Estimate of error standard deviation. If NULL (default), this is estimated from the mean squared residual
of the full fit (when n>p). If n le p, user should supply a value, using, for example the estimateSigma
function.
}
  \item{alpha}{
Significance level for  selection intervals (target is miscoverage alpha/2 in each tail)
}
  \item{verbose}{
Print out more details along the way.}
  \item{compute.si}{
Compute selection intervals? Default TRUE. Can slow down the computation, if only p-values are desired. 
}
 \item{tol.beta}{
Tolerance for determining if  coefficient is zero (default 1e-5)
}
  \item{tol.kkt}{
Tolerance for determining if  subgradient is zero (default 0.05)
}
  \item{one.sided}{
One-sided p-values? Default TRUE.}
}
\details{
This function computes selective p-values and confidence intervals (selection intervals) for the lasso
with a fixed value of the regularization parameter lambda.
The predictors and response variable are first centered. The selection interval estimation involves numerical search and can be fragile:
if the observed statistics are too close to either end of the truncation interavls (Vm, Vp) (see ref), the interval of desired coverage
cannot be computed and is set to (-Inf, +Inf). The output \code{tailarea} gives the achieved Gaussian tail area
for reported  interval--- these should be close to alpha/2.
}
\value{  
\item{pv}{p-values for each predictor}
\item{ci}{Selection intervals}
\item{tailarea}{Realized tail areas (lower and upper) for each  selection interval}
\item{eta}{Linear functionals defining each parameter}
\item{vm}{Lower bounds for functionals}
\item{vp}{Upper bounds for functionals}
\item{pred}{Predictors in active set}
\item{alpha}{Desired coverage (alpha/2 in each tail)}
\item{sigma}{Value of error standard deviation (sigma) used}
\item{one.sided}{Are the pvalues one-sided}
\item{lambda}{Value of regularization parameter lambda used}
\item{call}{The call to fixedLassoInf}
}

\references{
 Jason D. Lee, Dennis L. Sun, Yuekai Sun, Jonathan E. Taylor (2014). Exact post-selection inference, with application to the lasso. arXiv:1311.6238
}
\author{
Joshua Loftus, Stephen Reid, Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
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
# first run  glmnet
#gfit=glmnet(x,y,standardize=F)
#lam = .1
#extract coef for a given lam; Note the 1/n factor!
#bhat = coef(gfit, s=lam/n, exact=TRUE)[-1]

# compute fixed lambda p-values and selection intervals
#aa=fixedLassoInf(x,y,bhat,lam,sigma=sigma)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

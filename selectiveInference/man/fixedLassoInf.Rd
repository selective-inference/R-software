\name{fixedLassoInf}
\alias{fixedLassoInf}

\title{
Inference for the lasso, with a fixed lambda 
}
\description{
Compute p-values and confidence intervals for the lasso estimate, at a 
fixed value of the tuning parameter lambda
}
\usage{
fixedLassoInf(x, y, beta, lambda, intercept=TRUE, sigma=NULL, alpha=0.1,
         type=c("partial","full"), tol.beta=1e-5, tol.kkt=0.1,
         gridrange=c(-100,100), gridpts=10000, verbose=FALSE) 
}
\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
  \item{beta}{
Estimated lasso coefficients (e.g., from glmnet). This is of length p
(so the intercept is not included as the first component).

 Be careful! This function uses the "standard" lasso objective
  \deqn{
    1/2 \|y - x \beta\|_2^2 + \lambda \|\beta\|_1.
  }
 In contrast, glmnet multiplies the first term by a factor of 1/n.
 So after running glmnet, to extract the beta corresponding to a value lambda, 
 you need to use \code{beta = coef(obj, s=lambda/n)[-1]},
 where obj is the object returned by glmnet (and [-1] removes the intercept,
 which glmnet always puts in the first component)
}      
  \item{lambda}{
Value of lambda used to compute beta. See the above warning
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
\item{intercept}{
Was glmnet run with an intercept=TRUE? Default is TRUE
}
\item{type}{Contrast type for p-values and confidence intervals: default is
"partial"---meaning that the contrasts tested are the partial population 
regression coefficients, within the active set of predictors; the alternative is
"full"---meaning that the full population regression coefficients are tested.
The latter does not make sense when p > n.}
\item{tol.beta}{
Tolerance for determining if a coefficient is zero
}
\item{tol.kkt}{
Tolerance for determining if an entry of the subgradient is zero 
}
\item{gridrange}{
Grid range for constructing confidence intervals, on the standardized scale
}
\item{gridpts}{
Number of grid points for constructing confidence intervals, within the specified
grid range}
\item{verbose}{
Print out progress along the way? Default is FALSE}
}

\details{
This function computes selective p-values and confidence intervals for the lasso,
given a fixed value of the tuning parameter lambda. The confidence interval 
construction involves numerical search and can be fragile:
if the observed statistics are too close to either end of the truncation intervals 
(vlo and vup, see references), then the interval of desired coverage
cannot be computed and is set to [-Inf, +Inf]. The output \code{tailarea} gives 
the achieved Gaussian tail area for reported interval---these should be close to alpha/2,
and can be used for error-checking purposes. 

NOTE: typically it is best to standardize the features x before running lasso. This function  standardizes  so that colMeans(x)=0, colMeans(x^2)=1 (at least when intercept=TRUE).  We recommend  standardizing the data in this way before running the L1 fitter (eg glmnet) and  fixedLassoInf.
 You can do this via   x <- scale(x,T,T)/(sqrt(nrow(x))-1). See the examples below.
  If you choose not to scale the features,
you need to at last center them  via  x <- scale(x,T,F).
}
\value{  
\item{type}{Type of coefficients tested (partial or fall)}
\item{lambda}{Value of tuning parameter lambda used}
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
\item{call}{The call to lassoInf}
}

\references{
Jason Lee, Dennis Sun, Yuekai Sun, Jonathan Taylor (2014). 
Exact post-selection inference, with application to the lasso. arXiv:1311.6238.
}
\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Stephen Reid}

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
#lambda = .1
#extract coef for a given lambda; Note the 1/n factor!
#beta = coef(gfit, s=lambda/n, exact=TRUE)[-1]

# compute fixed lambda p-values and selection intervals
#aa=fixedLassoInf(x,y,beta,lambda,sigma=sigma)
#
# as above, but use lar function instead to get initial lasso fit (should get same result)
# fit=lar(x,y,normalize=F)
# beta=coef(fit,s=lambda,mode="lambda")
# fixedLassoInf(x,y,beta,lambda,sigma=sigma)


##
# states data example: predict life expectancy
#x=state.x77[,-4]
#y=state.x77[,4]
#x=scale(x,T,T)
#n=nrow(x)
#  estimate sigma from cross-validated lasso fit
#  cvf=cv.glmnet(x,y)
# sigmahat=estimateSigma(x,y,stand=F)$sigmahat
#
# fit lasso
# lambda=n*cvf$lambda.min
# gfit=glmnet(x,y,standardize=F)
# bhat=coef(gfit, s=lambda/n, exact=TRUE)[-1]
#
# compute p-values and confidence intervals
# fixedLassoInf(x,y,bhat,lambda,sigma=sigmahat)
}

\name{fixedLassoInf}
\alias{fixedLassoInf}

\title{
Inference for fixed lambda lasso
}
\description{
Compute p-values and selection intervals for lasso fit with fixed value of regularization parameter
lambda
}
\usage{
fixedLassoInf(x, y, bhat, lam, sigma, alpha = 0.1, trace = F, compute.si = TRUE, tol = 1e-05, one.sided = TRUE)
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
  \item{sigma}{
Estimate of error standard deviation
}
  \item{alpha}{
Significance level for  selection intervals (target is miscoverage alpha/2 in each tail)
}
  \item{trace}{
Print out more details along the way.}
  \item{compute.si}{
Compute selection intervals? Default TRUE. Can slow down the computation, if only p-values are desired. 
}
  \item{tol}{
Tolerance for constraints (default 1e-5)
}
  \item{one.sided}{
One-sided p-values? Default TRUE.}
}
\details{
This function computes selective p-values and confidence intervals (selection intervals) for the lasso
with a fixed value of the regularization parameter lambda. The selection interval estimation involves numerical search and can be fragile:
if the observed statistics are too close to either end of the truncation interavls (Vm, Vp) (see ref), the interval of desired coverage
cannot be computed and is set to (-Inf, +Inf). The output \code{tailarea} gives the achieved Gaussian tail area
for reported  interval--- these should be close to alpha/2.
}
\value{  
\item{pv}{p-values for each predictor}
\item{ci}{Selection intervals}
\item{tailarea}{Realized tail areas (left and right) for each  selection interval}
\item{eta}{Linear functionals defining each parameter}
\item{vm}{Lower bounds for functionals}
\item{vp}{Upper bounds for functionals}
\item{pred}{Predictors in active set}
\item{alpha}{Desired coverage (alpha/2 in each tail)}
\item{sigma}{Value of error standard deviation (sigma) used}
\item{one.sided}{Are the pvalues one-sided}
\item{lambda}{Value of regularization parameter lambda used}
\item{call}{The call to forwardStepInf}
}

\references{
 Jason D. Lee, Dennis L. Sun, Yuekai Sun, Jonathan E. Taylor (2014). Exact post-selection inference, with application to the lasso. arXiv:1311.6238
}
\author{
Joshua Loftus, Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

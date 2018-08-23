\name{randomizedLassoInf}
\alias{randomizedLassoInf}

\title{
Inference for the randomized lasso, with a fixed lambda 
}
\description{
Compute p-values and confidence intervals based on selecting
an active set with the randomized lasso, at a 
fixed value of the tuning parameter lambda and using Gaussian
randomization.
}
\usage{
randomizedLassoInf(rand_lasso_soln, 
		   targets=NULL,
                   level=0.9,
                   sampler=c("norejection", "adaptMCMC"),
                   nsample=10000,
                   burnin=2000,
                   opt_samples=NULL)
}
\arguments{
\item{rand_lasso_soln}{
A randomized lasso solution as returned by \code{randomizedLasso}.
}
\item{targets}{
If not NULL, should be a list with entries \code{observed_target, cov_target, crosscov_target_internal}.
The \code{observed_target} should be (pre-selection) asymptotically Gaussian around targeted
parameters. The quantity \code{cov_target} should be an estimate of the (pre-selection) covariance
of \code{observed_target}. Finally, \code{crosscov_target_internal} should be an estimate of
the (pre-selection) covariance of \code{observed_target} and the internal representation of the
data of the LASSO. For both \code{"gaussian"} and \code{"binomial"}, this is the vector
  \deqn{
\hat{\beta}_{E,MLE}, X_{-E}^T(y - \mu(X_E\hat{\beta}_{E,MLE}))
  } 
For example, this cross-covariance could be estimated by jointly bootstrapping the target 
of interest and the above vector.
}
\item{level}{
Level for confidence intervals.
}
\item{sampler}{
Which sampler to use -- default is a no-rejection sampler. Otherwise
use MCMC from the adaptMCMC package.
}
\item{nsample}{
Number of samples of optimization variables to sample.
}
\item{burnin}{
How many samples of optimization variable to discard (should be less than nsample).
}
\item{opt_samples}{
Optional sample of optimization variables. If not NULL then no MCMC will be run.
}
}

\details{
This function computes selective p-values and confidence intervals for a
randomized version of the lasso,
given a fixed value of the tuning parameter lambda. 

}
\value{  
\item{targets}{A list with entries \code{observed_target, cov_target, crosscov_target_internal}. See argument description above.}
\item{pvalues}{P-values testing hypotheses that each specific target is 0.}
\item{ci}{Confidence interval for parameters determined by \code{targets}.}
}

\references{
Jelena Markovic and Jonathan Taylor (2016).
Bootstrap inference after using multiple queries for model selection. arxiv.org:1612.07811

Xiaoying Tian and Jonathan Taylor (2015).
Selective inference with a randomized response. arxiv.org:1507.06739

Xiaoying Tian, Snigdha Panigrahi, Jelena Markovic, Nan Bi and Jonathan Taylor (2016).
Selective inference after solving a convex problem. 
arxiv.org:1609.05609

}
\author{Jelena Markovic, Jonathan Taylor}

\examples{
set.seed(43)
n = 50
p = 10
sigma = 0.2
lam = 0.5

X = matrix(rnorm(n*p), n, p)
X = scale(X, TRUE, TRUE) / sqrt(n-1)

beta = c(3,2,rep(0,p-2))
y = X\%*\%beta + sigma*rnorm(n)

result = randomizedLasso(X, y, lam)
inf_result = randomizedLassoInf(result)
}
 

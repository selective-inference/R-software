\name{fsInf_maxZ}
\alias{fsInf_maxZ}
\title{
Selective inference for forward stepwise regression
}
\description{
Computes maxZ selective p-values and confidence intervals for forward 
stepwise regression
}
\usage{
                    
fsInf_maxZ(obj, sigma=NULL, alpha=0.1, k=NULL,   ndraw=8000, burnin=2000,verbose=FALSE)
    
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
\item{ndraw}{Number of Monte Carlo samples generated}
\item{burnin}{
Number of samples discarded at the beginning of the chain
}
\item{verbose}{Print out progress along the way? Default is FALSE}
}

\details{
This function computes selective maxZ p-values
for forward stepwise regression. These p-values are independent the under null,
so that stopping via the forwardStop rule yields guaranteed FDR control
}

\value{  
\item{pv}{P-values for each model in the sequence}
\item{k}{Value of k specified in call}
\item{khat}{When type is "active", this is an estimated stopping point
declared by \code{\link{forwardStop}}}
\item{sigma}{Value of error standard deviation (sigma) used}
\item{vars}{Variables in active set}
\item{sign}{Signs of active coefficients}
\item{alpha}{Desired coverage (alpha/2 in each tail)}
\item{realized_maxZ}{Value of maxZ statistic computed at each step}
\item{call}{The call to fsInf_maxZ}
}

\references{
Will Fithian, Jonathan Taylor,  Ryan Tibshirani, and Rob Tibshirani (2015).
Selective sequential model selection.  arXiv:1512.02565..


}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid}

\seealso{\code{\link{fs}}}

\examples{
set.seed(33)
n = 50
p = 10
sigma = 1
x = matrix(rnorm(n*p),n,p)
beta = c(3,2,rep(0,p-2))
y = x\%*\%beta + sigma*rnorm(n)

# run forward stepwise
fsfit = fs(x,y)

# compute sequential p-values a
# (sigma estimated from full model)
out.seq = fsInf_maxZ(fsfit)
out.seq
}

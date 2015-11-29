\name{manyMeans}
\alias{manyMeans}
\title{
Selective inference for many normal means
}
\description{
Computes p-values and confidence intervals for the largest k
among many normal means
}
\usage{
manyMeans(y, alpha=0.1, bh.q=NULL, k=NULL, sigma=1, verbose=FALSE)
}

\arguments{
  \item{y}{Vector of outcomes (length n)}
\item{alpha}{
Significance level for confidence intervals (target is miscoverage alpha/2 in each tail)
}
\item{bh.q}{q parameter for BH(q) procedure}
\item{k}{Number of means to consider}
\item{sigma}{Estimate of error standard deviation}
\item{verbose}{Print out progress along the way? Default is FALSE}
}
\details{
This function compute p-values and confidence intervals for the largest k 
among many normal means. One can specify a fixed number of means k to consider, 
or choose the number to consider via the BH rule.
}
  
\value{  
\item{mu.hat}{ Vector of length n containing the estimated signal sizes.
If a sample element is not selected, then its signal size estimate is 0}
\item{selected.set}{Indices of the  vector y of the sample elements that 
were selected by the procedure (either BH(q) or top-K). Labelled "Selind" in output table.}
\item{pv}{P-values for selected signals}
\item{ci}{Confidence intervals}
\item{method}{Method used to choose number of means}
\item{sigma}{Value of error standard deviation (sigma) used}
\item{bh.q}{BH q-value used}
\item{k}{Desired number of means}
\item{threshold}{Computed cutoff}
\item{call}{The call to manyMeans}
}

\references{
Stephen Reid, Jonathan Taylor, and Rob Tibshirani (2014).
Post-selection point and interval estimation of signal sizes in Gaussian samples.
arXiv:1405.3340.
}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid}

\examples{
set.seed(12345)
n = 100 
mu = c(rep(3,floor(n/5)), rep(0,n-floor(n/5))) 
y = mu + rnorm(n)
out = manyMeans(y, bh.q=0.1)
out
}

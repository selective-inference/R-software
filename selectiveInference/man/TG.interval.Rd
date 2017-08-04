\name{TG.interval}
\alias{TG.interval}

\title{
Truncated Gaussian confidence interval.
}
\description{
Compute truncated Gaussian interval of Lee et al. (2016) with
arbitrary affine selection and covariance.
Z should satisfy A %*% Z elementwise less then or equal b.
}
\usage{
TG.interval(Z, A, b, eta, Sigma=NULL, alpha=0.1, 
	   gridrange=c(-100,100),
           gridpts=100, 
           griddepth=2, 
           flip=FALSE, 
           bits=NULL) 
}
\arguments{
\item{Z}{
Observed data (assumed to follow N(mu, Sigma) with sum(eta*mu)=null_value)
}
\item{A}{
Matrix specifiying affine inequalities AZ <= b
}
\item{b}{
Offsets in the affine inequalities AZ <= b.
}      
\item{eta}{
Determines the target sum(eta*mu) and estimate sum(eta*Z).
} 
\item{Sigma}{
Covariance matrix of Z. Defaults to identity.
} 
\item{alpha}{
Significance level for confidence intervals (target is miscoverage alpha/2 in each tail)
}
\item{gridrange}{
Grid range for constructing confidence intervals, on the standardized scale.
}
\item{gridpts}{
???????
}
\item{griddepth}{
???????
}
\item{flip}{
???????
}
\item{bits}{
Number of bits to be used for p-value and confidence interval calculations. Default is 
NULL, in which case standard floating point calculations are performed. When not NULL, 
multiple precision floating point calculations are performed with the specified number 
of bits, using the R package \code{Rmpfr} (if this package is not installed, then a 
warning is thrown, and standard floating point calculations are pursued).
Note: standard double precision uses 53 bits
so, e.g., a choice of 200 bits uses about 4 times double precision. The confidence
interval computation is sometimes numerically challenging, and the extra precision can be 
helpful (though computationally more costly). In particular, extra precision might be tried 
if the values in the output columns of \code{tailarea} differ noticeably from alpha/2.
}
}

\details{
This function computes selective confidence intervals based on the polyhedral
lemma of Lee et al. (2016). 
}

\value{  
\item{int}{Selective confidence interval.}
\item{tailarea}{Realized tail areas (lower and upper) for each confidence interval.}
}

\references{
Jason Lee, Dennis Sun, Yuekai Sun, and Jonathan Taylor (2016). 
Exact post-selection inference, with application to the lasso. Annals of Statistics, 44(3), 907-927.

Jonathan Taylor and Robert Tibshirani (2017) Post-selection inference for math L1-penalized likelihood models.
Canadian Journal of Statistics, xx, 1-21. (Volume still not posted)

}
\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid}

\examples{
A = diag(5)
b = rep(1, 5)
Z = rep(0, 5)
Sigma = diag(5)
eta = as.numeric(c(1, 1, 0, 0, 0))
TG.interval(Z, A, b, eta, Sigma)
}
 
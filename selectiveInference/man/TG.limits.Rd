\name{TG.limits}
\alias{TG.limits}

\title{
Truncation limits and standard deviation.
}
\description{
Compute truncated limits and SD for use in computing
p-values or confidence intervals of Lee et al. (2016).
Z should satisfy A %*% Z elementwise less then or equal b.
}
\usage{
TG.limits(Z, A, b, eta, Sigma)  
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
}
\details{
This function computes the limits of truncation and the implied
standard deviation in the polyhedral lemma of Lee et al. (2016). 
}

\value{  
\item{vlo}{Lower truncation limits for statistic}
\item{vup}{Upper truncation limits for statistic}
\item{sd}{Standard error of sum(eta*Z)}
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
TG.limits(Z, A, b, eta, Sigma)
}
 
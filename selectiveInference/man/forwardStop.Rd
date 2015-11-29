\name{forwardStop}
\alias{forwardStop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ForwardStop rule for sequential p-values
}
\description{
Computes the ForwardStop sequential stopping rule of G'Sell et al (2014)
}
\usage{
forwardStop(pv, alpha=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pv}{
Vector of **sequential** p-values, for example from fsInf or larInf
}
  \item{alpha}{
Desired type FDR level (between 0 and 1)
}
}
\details{
Computes the ForwardStop sequential stopping rule of G'Sell et al (2014).
Guarantees FDR control at the level alpha, for independent p-values.
}
\value{
Step number for sequential stop.
}
\references{
Max Grazier G'Sell, Stefan Wager, Alexandra Chouldechova, and Rob Tibshirani (2014).
Sequential selection procedures and Fflse Discovery Rate Control. arXiv:1309.5352.
To appear in Journal of the Royal Statistical Society: Series B.
}
\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid}

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

# compute sequential p-values and confidence intervals
# (sigma estimated from full model)
out = fsInf(fsfit)
out

# estimate optimal stopping point
forwardStop(out$pv, alpha=.10)
}

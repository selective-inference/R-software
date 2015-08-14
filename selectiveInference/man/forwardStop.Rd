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
Max Grazier G'Sell, Stefan Wager, Alexandra Chouldechova, Rob Tibshirani (2014).
Sequential selection procedures and Fflse Discovery Rate Control. arXiv:1309.5352.
To appear in Journal of the Royal Statistical Society: Series B.
}
\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Stephen Reid}

\examples{
#NOT RUN
#set.seed(433)
#n=50
#p=10
#sigma=.7
#x=matrix(rnorm(n*p),n,p)

#x=scale(x,T,T)/sqrt(n-1)
#beta=c(4,2,0,0,rep(0,p-4))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#first run forward stepwise
# fsfit=fs(x,y)
#
# aa=fsInf(fsfit)
# forwardStop(aa, alpha=.10)
}

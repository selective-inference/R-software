\name{forwardStop}
\alias{forwardStop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forward stop rule for sequential inference
}
\description{
Computes the forward Stop sequential stopping rule  of  G'Sell et al (2014).
}
\usage{
forwardStop(pv, alpha = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pv}{
vector of **sequential** p-values, for example from forwardStepInf or larInf
}
  \item{alpha}{
Desired type FDR level (between 0 and 1)
}
}
\details{
Computes the forward Stop sequential stopping rule  of  G'Sell et al (2014).
Guarantees FDR control at level alpha.
}
\value{
Step number for sequential stop.
}
\references{
Max Grazier G'Sell, Stefan Wager, Alexandra Chouldechova, Robert Tibshirani (2014).
Sequential Selection Procedures and False Discovery Rate Control. arXiv:1309.5352.
To appear JRSSB.
}
\author{
Joshua Loftus, Stephen Reid,  Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
}


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (pv, alpha = 0.1) 
{
    val = -(1/(1:length(pv))) * cumsum(log(1 - pv))
    oo = which(val <= alpha)
    if (length(oo) == 0) 
        out = 0
    if (length(oo) > 0) 
        out = oo[length(oo)]
    return(out)
  }
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

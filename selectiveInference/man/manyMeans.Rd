\name{manyMeans}
\alias{manyMeans}

\title{
P-values and selection intervals  for the  largest k
among many normal means
}
\description{
This function compute p-values and selection intervals  for the  largest k
among many normal means.
}
\usage{
manyMeans(y,  alpha = 0.05, bh.q=NULL, k=NULL,sigma=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Vector of outcomes (length n)}
 \item{alpha}{
Significance level for  selection intervals (target is miscoverage alpha/2 in each tail)
}
\item{bh.q}{q parameter for BH(q) procedure. Default: NULL}
\item{k}{Number of largest elements to consider. Default: NULL}
  \item{sigma}{Estimate of error standard deviation}
}
\details{
This function compute p-values and selection intervals  for the  largest k 
among many normal means. One can specify a fixed number of means or  choose the number by the BH rule.
}

  
\value{  
\item{muhat}{ Vector of length length(y) containing the estimated signal size. If a sample element is not selected, then its signal size estimate is 0}
\item{selected set}{Indices into the vector y of the sample elements that were selected by our procedure (either BH(q) or top-K)}
\item{selint}{Matrix with two columns and number of rows equal to number of elements in selected.set. Provides the post-selection  bounds for the estimated signal sizes of selected elements. Selection intervals given in rows in the same order as encountered in selected.set}
\item{p.vals}{p-values for each predictor}
\item{method}{Method used to choose number of means}
\item{q}{BH q-value used}
\item{k}{Desired number of means}
\item{threshold}{Computed cutoff}
\item{call}{The call to manyMeans}
}

\references{
Stephen Reid, Jonathan Taylor, Robert Tibshirani (2014).
Post-selection point and interval estimation of signal sizes in Gaussian samples.
arXiv:1405.3340
}
\author{
Joshua Loftus, Stephen Reid, Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{


#NOT RUN
#set.seed(12345)
#n = 100 # sample size
#signal = 3 # signal size
#mu = c(rep(signal, floor (n/5)), rep(0, n-floor(n/5))) # 20% of elements get the signal; rest 0
#y = mu + rnorm (n, 0, 1)
#mmObj = manyMeans(y, bh.q=0.1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

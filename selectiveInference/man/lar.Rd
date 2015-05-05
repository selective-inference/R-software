\name{lar}
\alias{lar}
\title{
Least angle regression
}
\description{
Computes a least angle regression
}
\usage{
lar(x, y, maxsteps = 2000, minlam = 0, verbose = FALSE, 
    normalize = TRUE, intercept = TRUE, eps = .Machine$double.eps) 
}

\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
  \item{maxsteps}{
Maximum number of steps to take (default=min(n,p,2000))
}
  \item{minlam}{
 Smallest lambda value to consider
}
\item{verbose}{Print out details along the way? Default: FALSE}
\item{normalize}{Center x and y, and standardize columns of x to have variance 1? Default TRUE}
\item{intercept}{Include an intercept? Default TRUE}
\item{eps}{Numerical precision. Default: .Machine$double.eps}
}
\details{
LARS is described in detail in Efron, Hastie, Johnstone and
     Tibshirani (2002). This function is very similar to the lars function in the
 LARS library in R, but
returns additional information for the selective inference computations.
}




\value{
\item{lambda}{Values of lambda considered.}
\item{action}{Predictors entered in order}
\item{df}{Degrees of freedom at each step}
\item{beta}{Matrix of estimated coefficients}
    \item{completepath}{}
\item{Gamma}{Constraint matrix for inference computation}
\item{nk}{Number of constraints for inference at each knot}
 \item{mp}{Value of M+ (for internal use)}
\item{mu}{Estimated intercept}
\item{meanx}{Means of each column of x}
\item{normx}{Norms of each column of x}
\item{normalize}{Was data normalized?}
\item{intercept}{Is an intercept included in model? }
\item{type}{Call type: currently only "LAR" allowed}
\item{call}{The call to forwardStep}
}

\references{
Efron, Hastie, Johnstone and Tibshirani (2002) "Least Angle
     Regression" (with discussion) _Annals of Statistics_; see also
     <URL:
     http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf>.

     Hastie, Tibshirani and Friedman (2002) Elements of Statistical
     Learning, Springer, NY.
}
\author{
Joshua Loftus, Stephen Reid, Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{lar}}
}
\examples{
#NOT RUN
#set.seed(33)
#n=20
#p=10
#sigma=1
#x=matrix(rnorm(n*p),n,p)
#generate data
#beta=c(3,3,rep(0,p-2))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#run lar
#larfit=lar(x,y,verbose=TRUE)
                                      
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

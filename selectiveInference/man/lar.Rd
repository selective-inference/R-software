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
     intercept = TRUE, normalize = TRUE) 
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
\item{intercept}{Include an intercept? Default TRUE}
\item{normalize}{Center x and y, and standardize columns of x to have variance 1? Default TRUE}
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
\item{sign}{Sign of each coefficient as predictor is entered}
\item{df}{Degrees of freedom at each step}
\item{beta}{Matrix of estimated coefficients}
    \item{completepath}{Was the complete LAR path computed?}
\item{bls}{LAR coefficients for the last model}
\item{Gamma}{Constraint matrix for inference computation}
\item{nk}{Number of constraints for inference at each knot}
 \item{mp}{Value of M+ (for internal use)}
\item{bx}{Means of each column of x}
\item{by}{Mean of y}
\item{sx}{Standard deviation of each column of x}
\item{intercept}{Is an intercept included in model? }
\item{normalize}{Was data normalized?}
\item{call}{The call to forwardStep}
}

\references{
Efron, Hastie, Johnstone and Tibshirani (2002) "Least Angle
     Regression" (with discussion) _Annals of Statistics_; see also
     <URL:
     http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf>.

     Hastie, Tibshirani and Friedman (2002, 2009) Elements of Statistical
     Learning, Springer, NY.
}

\author{Ryan Tibshirani, Rob Tibshirani,  Jonathan Taylor, Max G'Sell, Stephen Reid}


\seealso{
 \code{\link{predict.lar}}
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
#plot(larfit)                                      
}

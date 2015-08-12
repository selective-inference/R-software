\name{predict.fs}
\alias{predict.fs}
\title{
Prediction from forward stepwise fit
}
\description{
Make predictions or extract coefficients from a fitted forward stepwise model
}
\usage{
\method{predict}{fs}(object, newx, s, ...)
}

\arguments{
  \item{object}{
Object returned by a call to \code{fs} function
}
  \item{newx}{
          The matrix of x values at which the
          fit is required. 
}
  \item{s}{
        a step number, or vector of step values, indexing the path.
}
\item{\dots}{Any additional arguments for fs}
}
  
\value{
Either a vector/matrix of fitted values, or a vector/matrix of
     coefficients.
}
\references{

     Hastie, Tibshirani and Friedman (2002) Elements of Statistical
     Learning, Springer, NY.

      Jonathan Taylor, Richard Lockhart,  Ryan Tibshirani and Rob Tibshirani (2015) 
    Exact Post-selection Inference for Forward Stepwise and Least Angle Regression.
arXiv:1401.3889
}

\author{Ryan Tibshirani, Rob Tibshirani,  Jonathan Taylor, Max G'Sell, Stephen Reid}


\examples{
#NOT RUN
#set.seed(33)
#n=200
#p=20
#sigma=1
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,T)/sqrt(n-1)
#generate y
#beta=c(3,-2,rep(0,p-2))
#beta=c(rep(3,10),rep(0,p-10))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#obj=fs(x,y)
#fit=predict.fs(obj,x)
}

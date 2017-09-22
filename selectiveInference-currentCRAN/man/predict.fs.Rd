\name{predict.fs}
\alias{predict.fs}
\alias{coef.fs}

\title{
Prediction and coefficient functions for forward stepwise
regression
}
\description{
Make predictions or extract coefficients from a forward stepwise object
}
\usage{
\method{predict}{fs}(object, newx, s, ...)
\method{coef}{fs}(object, s, ...)
}

\arguments{
\item{object}{
Object returned by a call to \code{fs} function
}
\item{newx}{
Matrix of x values at which the predictions are desired. If NULL,
the x values from forward stepwise fitting are used
}
\item{s}{
Step number(s) at which predictions or coefficients are desired
}
\item{\dots}{Additional arguments}
}
  
\value{
Either a vector/matrix of predictions, or a vector/matrix of coefficients.
}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid}

\examples{
set.seed(33)
n = 200
p = 20
sigma = 1
x = matrix(rnorm(n*p),n,p)
beta = c(rep(3,10),rep(0,p-10))
y = x\%*\%beta + sigma*rnorm(n)

# run forward stepwise and predict functions
obj = fs(x,y)
fit = predict(obj,x,s=3)
}

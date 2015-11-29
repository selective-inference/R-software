\name{predict.lar}
\alias{predict.lar}
\alias{coef.lar}

\title{
Prediction and coefficient functions for least angle regression
}
\description{
Make predictions or extract coefficients from a least angle regression object
}
\usage{
\method{predict}{lar}(object, newx, s, mode=c("step","lambda"), ...)
\method{coef}{lar}(object, s, mode=c("step","lambda"), ...)
}

\arguments{
\item{object}{
Object returned by a call to \code{lar} function
(not the \code{lars} function!)
}
\item{newx}{
Matrix of x values at which the predictions are desired. If NULL,
the x values from least angle regression fitting are used
}
\item{s}{
Step number(s) or lambda value(s) at which predictions or coefficients 
are desired
}
\item{mode}{Either "step" or "lambda", determining the role of s (above)}

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

# run lar and predict functions
obj = lar(x,y)
fit = predict(obj,x,s=3)
}

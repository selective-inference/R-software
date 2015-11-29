\name{estimateSigma}
\alias{estimateSigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the noise standard deviation in regression 
}
\description{
Estimates the standard deviation of the noise, for use in the selectiveInference 
package
}
\usage{
estimateSigma(x, y, intercept=TRUE, standardize=TRUE)
}
\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
\item{intercept}{Should glmnet be run with an intercept? Default is TRUE}
\item{standardize}{Should glmnet be run with standardized predictors? Default is TRUE}
}
\details{
This function estimates the standard deviation of the noise, in a linear regresion setting.
A lasso regression is fit, using cross-validation to estimate the tuning parameter lambda.
With sample size n, yhat equal to the predicted values and df being the number of nonzero 
coefficients from the lasso fit, the estimate of sigma is \code{sqrt(sum((y-yhat)^2) / (n-df-1))}.
Important: if you are using glmnet to compute the lasso estimate, be sure to use the settings
for the "intercept" and "standardize" arguments in glmnet and estimateSigma. Same applies to fs 
or lar, where the argument for standardization is called "normalize".
}
\value{
\item{sigmahat}{The estimate of sigma}
\item{df}{The degrees of freedom of lasso fit used}
}
\references{
Stephen Reid, Jerome Friedman, and Rob Tibshirani (2014). 
A study of error variance estimation in lasso regression. arXiv:1311.5274.
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

# estimate sigma
sigmahat = estimateSigma(x,y)$sigmahat

# run sequential inference with estimated sigma
out = fsInf(fsfit,sigma=sigmahat)
out
}



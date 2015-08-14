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
Important: if you are using glmnet to compute the lasso estimate, be sure to use the same values
for the intercept and standardized argument in glmnet and estimateSigma. Same applies to  fs or lar,
where the argument for standardization is called "normalize"
}
\value{
\item{sigmahat}{The estimate of sigma}
\item{df}{The degrees of freedom of lasso fit used}
}
\references{
Stephen Reid, Jerome Friedman and Rob Tibshirani (2014). 
A study of error variance estimation in lasso regression. arXiv:1311.5274.
}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Stephen Reid}

\examples{
#NOT RUN
#set.seed(43)
#n=50
#p=10
#sigma=.7
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,F)
#beta=c(3,2,0,0,rep(0,p-4))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)

#out=estimateSigma(x,y)
#a=fs(x,y)
#out2=fsInf(a,sigma=out$sigmahat)

}



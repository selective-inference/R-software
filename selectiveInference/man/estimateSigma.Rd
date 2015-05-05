\name{estimateSigma}
\alias{estimateSigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forward stepwise regression
}
\description{
This function estimates the standard deviaiton of the noise, for use in the selectiveInference package
}
\usage{
estimateSigma(x, y)
}
\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
}
\details{
This function estimates the standard deviation of the noise, in a linear regresion setting.
A lasso regression is fit, using cross-validation to estimate the tuning parameter lambda.
With sample size n,  yhat equal to the predicted values and df being the number of non-zero coefficients from
the lasso fit, the estimate of sigma is  sqrt( sum( (y-yhat)^2)/(n-df-1))
}
\value{
\item{sigmahat}{The estimate of sigma}
\item{df}{The estimate degrees of freedom of lasso fit ued for estimating sigma}
}
\references{
Stephen Reid, Jerome Friedman and  Rob Tibshirani (2014). A Study of Error Variance Estimation in Lasso Regression. arXiv:1311.5274.
}

\author{
Joshua Loftus, Stephen Reid, Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
}


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
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

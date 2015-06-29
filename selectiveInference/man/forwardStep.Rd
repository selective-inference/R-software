\name{forwardStep}
\alias{forwardStep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forward stepwise regression
}
\description{
This function implements forward stepwise regression, for use in the selectiveInference package
}
\usage{
forwardStep(x, y, sigma=NULL, intercept=TRUE, nsteps = min(nrow(x)-1, ncol(x)),trace=FALSE)
}
\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
\item{sigma}{Optional estimate of standard deviation of the noise, for AIC.
 If NULL (default), this is estimated from the mean squared residual
of the full fit (when n>p). If n le p, user should supply a value, using, for example the estimateSigma
function.} 
\item{intercept}{Should an intercept be included on the model? Default TRUE}
  \item{nsteps}{
Number of steps to take (default=min(n,p))
}
\item{trace}{Printout steps along the way? Default: FALSE}
}
\details{
This function implements forward stepwise regression, adding the predictor that  maximizes the scaled inner product between the residual r and the 
orthogonalized predictors at each step. That is, argmax t(xr) %*% r/sqrt(diag(t(xr) %*% xr)),
where xr is the matrix of predictors with the current model removed.
If intercept is TRUE, the data x and y are first centered by the function.
}
\value{
\item{pred}{Predictors entered in order}
\item{scor}{Scaled inner products of entered predictors}
\item{bhat}{Partial regression coefficients of predictors entered at each step}
\item{rss}{Residual sum of squares after each predictor is entered}
\item{sigma}{Estimate of noise standard deviation, used for AIC}
\item{intercept}{The value of argument "intercept" specified in the call}
\item{aic}{AIC criterion after each predictor is added}
\item{aichat}{Step number minizing AIC (starting from the null model; first local minimum)}
\item{call}{The call to forwardStep}
}


\author{
Joshua Loftus, Stephen Reid, Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
}

\seealso{forwardStepInf}
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

#fsfit=forwardStep(x,y)
#
#out=forwardStepInf(fsfit,x,y)
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

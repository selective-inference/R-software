\name{fs}
\alias{fs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forward stepwise regression
}
\description{
This function implements forward stepwise regression, for use in the selectiveInference package
}
\usage{
fs(x, y, maxsteps=2000, verbose=FALSE, intercept=TRUE,
               normalize=TRUE)
}
\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}

  \item{maxsteps}{
Maximum number of steps to take (default=2000)
}
\item{verbose}{Print details along the way. Default: FALSE}

\item{intercept}{Should an intercept be included on the model? Default TRUE}
  \item{normalize}{Should features be standardized? Default TRUE}
}

\details{
This function implements forward stepwise regression, adding the predictor that  maximizes the scaled inner product between the residual r and the 
orthogonalized predictors at each step. That is, argmax t(xr) %*% r/sqrt(diag(t(xr) %*% xr)),
where xr is the matrix of predictors with the current model removed.
If intercept is TRUE, the data x and y are first centered by the function.
Note that fs uses the standard criterion for entering predictors, while  the \code{lars} function
with type=forward.stepwise uses a different, non-standard choice (correlation of residual with orthogonalized predictors)
}
\value{
\item{action}{Vector of predictors in order of entry}
\item{sign}{Signs of coefficients of entered predictors}
\item{df}{Degrees of freedom of each model}
\item{beta}{Matrix of regression coefficients for each model along the path; one model per column}
\item{completePath}{Was complete FS path computed?}
\item{bls}{Regression coefficients for the last model}
\item{Gamma}{Matrix that captures the polyhedral selection at each step}
\item{nk}{Number of constraints for each step in path}
\item{x}{Matrix of predictors used}
\item{y}{Vector of outcomes used}
\item{bx}{Vector of column means of x}
\item{by}{Mean of y}
\item{sx}{Standard deviation of each column of x}
\item{intercept}{The value of argument "intercept" specified in the call}
\item{normalize}{Was the data normalized?}
\item{call}{The call to fs}
}



\author{Ryan Tibshirani, Rob Tibshirani,  Jonathan Taylor, Max G'Sell, Stephen Reid}


\seealso{fsInf}
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

#fsfit=fs(x,y)
#
#out=fsInf(x,y)
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

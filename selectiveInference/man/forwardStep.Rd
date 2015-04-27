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
forwardStep(x, y, nsteps = min(nrow(x), ncol(x)))
}
\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
  \item{nsteps}{
Number of steps to take (default=min(n,p))
}
}
\details{
This function implements forward stepwise regression, adding the predictotr maximize the scaled inner product between the residual r and the 
orthogonalized predictors at each step. That is, argmax t(xr) %*% r/sqrt(diag(t(xr) %*% xr)).
The data x and y are first centered by the function.
}
\value{
\item{pred}{Predictors entered in order}
\item{scor}{Scaled inner products of entered predictors}
\item{bhat}{Partial regression coeficients of entered predictors}
\item{rss}{Residual sum of squares after each predictor is entered}
\item{sigma}{Estimate of noise standard deviation, used for AIC}
\item{aic}{AIC criterion afetr each predictor is added}
\item{aichat}{Step number minizing AIC (starting from the null model; first local minimum)}
\item{call}{The call to forwardStep}
}
\references{
%%
}

\author{
Joshua Loftus, Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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

#fsfit=forwardStep(x,y)
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

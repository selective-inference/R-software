\name{fs}
\alias{fs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forward stepwise regression
}
\description{
This function implements forward stepwise regression, for use in the 
selectiveInference package
}
\usage{
fs(x, y, maxsteps=2000, intercept=TRUE, normalize=TRUE, verbose=FALSE)
}
\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
  \item{maxsteps}{
Maximum number of steps to take 
}
\item{intercept}{Should an intercept be included on the model? Default is TRUE}
\item{normalize}{Should the predictors be normalized? Default is TRUE}
\item{verbose}{Print out progress along the way? Default is FALSE}
}

\details{
This function implements forward stepwise regression, adding the predictor that  
achieves the largest drop in RSS when added to the active model, at each step.
Note that this version of stepwise uses the more "standard" criterion for entering
variables, while that in the \code{lars} package (with type="step") uses a different
choice (entering the predictor that has the largest absolute correlation with the 
the residual).
}
\value{
\item{action}{Vector of predictors in order of entry}
\item{sign}{Signs of coefficients of predictors, upon entry}
\item{df}{Degrees of freedom of each active model}
\item{beta}{Matrix of regression coefficients for each model along the path; one model per column}
\item{completepath}{Was the complete stepwise path computed?}
\item{bls}{If completepath is TRUE, the full least squares coefficients}
\item{Gamma}{Matrix that captures the polyhedral selection at each step}
\item{nk}{Number of polyhedral constraints, at each step in path}
\item{x}{Matrix of predictors used}
\item{y}{Vector of outcomes used}
\item{bx}{Vector of column means of original x}
\item{by}{Mean of original y}
\item{sx}{Norm of each column of original x}
\item{intercept}{Was an intercept included?}
\item{normalize}{Were the predictors normalized?}
\item{call}{The call to fs}
}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Stephen Reid}

\seealso{
 \code{\link{fsInf}}, \code{\link{predict.fs}},\code{\link{coef.fs}}, \code{\link{plot.fs}}
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

#fsfit=fs(x,y)
#
#out=fsInf(x,y)
#
# states data example: predict life expectancy
#x=state.x77[,-4]
#y=state.x77[,4]
#x=scale(x,T,T)
#n=nrow(x)
#
#  estimate sigma from cross-validated lasso fit
#cvf=cv.glmnet(x,y)
#sigmahat=estimateSigma(x,y,stand=F)$sigmahat
#
# fit forwardStepwise
# fsfit=fs(x,y,normalize=F)
#
# compute p-values and confidence intervals
#  fsInf(fsfit,sigma=sigmahat)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

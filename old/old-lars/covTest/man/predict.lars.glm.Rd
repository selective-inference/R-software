\name{predict.lars.glm}
\alias{predict.lars.glm}
\title{
Function to make predictions from lars.glm fit
}
\description{
Function to make predictions from lars.glm fit
}
\usage{
\method{predict}{lars.glm}(object, x, lambda, time=NULL, status = NULL, ...)
}
\arguments{
  \item{object}{
}
Result of call to lars.glm
  \item{x}{
N by p matrix of predictors
}
  \item{lambda}{
Value of L1- regularization parameter at which predictions are desired
}
  \item{time}{
Optional N-vector of survival times, required for Cox Proportional hazards model.
}
  \item{status}{
Optional N-vector of censoring indicators, required for Cox Proportional hazards model. 1=failed; 0=censored.
}
 \item{...}{
       additional arguments (not used)
         }
}
\value{
Vector of predicted values, on the linear predictor scale.
}
\references{
Park, M.Y. and Hastie, T. (2007) 1l regularization path algorithm for generalized linear
models. JRSSB B 69(4), 659-677
}

\author{
Rob Tibshirani
}

\seealso{
lars.glm, covTest
}
\examples{
#logistic
x=matrix(rnorm(100*10),ncol=10)
x=scale(x,TRUE,TRUE)/sqrt(99)
y=4*x[,2]+rnorm(100)
y=1*(y>0)
a=lars.glm(x,y,family="binomial")
yhat=predict.lars.glm(a,x,family="binomial")

# Cox model
#y=6*x[,2]+rnorm(100)+10
#status=sample(c(0,1),size=length(y),replace=TRUE)
#a=lars.glm(x,y,status=status,family="cox")
#yhat=predict.lars.glm(a,x,family="cox")
}
\keyword{LARS}
\keyword{lasso}% __ONLY ONE__ keyword per line
\keyword{logistic}% __ONLY ONE__ keyword per line
\keyword{Cox model}% __ONLY ONE__ keyword per line

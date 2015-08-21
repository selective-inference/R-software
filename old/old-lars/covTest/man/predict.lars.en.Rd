\name{predict.lars.en}
\alias{predict.lars.en}
\title{
Function to make predictions from lars.en fit
}
\description{
Function to make predictions from lars.en fit
}
\usage{
\method{predict}{lars.en}(object, x, lambda, ...)
}
\arguments{
  \item{object}{
Result of call to lars.en
     }
  \item{x}{
N by p matrix of predictors
}
  \item{lambda}{
Value of L1- regularization parameter at which predictions are desired
}
  \item{...}{
       additional arguments (not used)
         }
}
\details{
Makes predictions from an objected returned by lars.en
}
\value{
Vector of predicted values.
}
\references{
Hui Zou and Trevor Hastie. Regularization and Variable Selection via the Elastic Net. JRSSB (2005) 67(2) 301-320. 
}

\author{
Rob Tibshirani
}


\seealso{
lars.en, covTest
}
\examples{
x=matrix(rnorm(100*10),ncol=10)
x=scale(x,TRUE,TRUE)/sqrt(99)
y=4*x[,2]+rnorm(100)
a=lars.en(x,y,lambda2=1)
yhat=predict.lars.en(a,x,.5)
}
\keyword{ LARS }
\keyword{ lasso }% __ONLY ONE__ keyword per line
\keyword{ elastic net }% __ONLY ONE__ keyword per line

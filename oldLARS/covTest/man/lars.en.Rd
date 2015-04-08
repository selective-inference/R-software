\name{lars.en}
\alias{lars.en}
\title{
Function to fit least angle regression path of solution for the elastic net.
}
\description{
Function to fit least angle regression path of solution for the elastic net.
}
\usage{
lars.en(x, y, lambda2,normalize=TRUE)
}
\arguments{
 \item{x}{
N by p matrix of predictors
}
  \item{y}{
N-vector of outcome values
}
  \item{lambda2}{
     Value of L2 penalty parameter 
}
  \item{normalize}{
Should columns of x be standardized?
}
}
\details{
This function estimates the  least angle regression path of solution for Ll-penalized (lasso) logistic regression
and the Cox proportional hazards model, using the R functions \code{enpath} and \code{coxpath}.
These latter functions use the predictor-corrector strategy devised by Park and Hastie (2007).
}
\value{
\item{beta}{Matrix whose rows of contain the estimated coefficients for each lambda value}
  \item{larsobj}{Result of call to lars on augmented data}
\item{mx}{Column means of x}
\item{sdx}{Column standard deviations of x}
\item{normalize}{Value of normalize argument in call to lars.en}
 \item{lambda}{Values of lambda used}
 \item{lambda2}{Value of lambda2 used}
  \item{act}{Actions (predictor added) at each step}
  \item{maxp}{Maximum number of predictors entered}
\item{call}{Call to lars.en}
}
\references{
Zou, H. and Hastie, Trevor (2005)
Regularization and Variable Selection via the Elastic Net. 
JRSSB 301-320,

Park, M. Y. & Hastie, T. (2007).
 l1-regularization path algorithm for generalized linear
models, Journal of the Royal Statistical Society Series B 69(4), 
}
\author{
Rob Tibshirani
}


\seealso{
predict.lars, covTest
}
\examples{
set.seed(432)
x=matrix(rnorm(100*10),ncol=10)
x=scale(x,TRUE,TRUE)/sqrt(99)
beta=c(3,rep(0,9))
y=x\%*\%beta+.4*rnorm(100)
a=lars.en(x,y,lambda2=.5)
  }
\keyword{ LARS }
\keyword{ lasso }% __ONLY ONE__ keyword per line
\keyword{ elastic net }% __ONLY ONE__ keyword per line

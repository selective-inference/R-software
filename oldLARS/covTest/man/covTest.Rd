\name{covTest}
\alias{covTest}
\title{
Compute the covariance test for adaptive linear modelling
}
\description{
This function  computes covariance test for inference in  adaptive linear modelling, for lasso (least angle
regression) in the Gaussian case, binomial/logistic and Cox proportional hazards survival models.
This package should be considered EXPERIMENTAL. The background paper is not yet published
and rigorous theory does not yet exist for the logistic and Cox models.
We have currently disabled the Cox option, as it is not yet reliable.

}
\usage{
covTest(fitobj, x, y, sigma.est = "full", status = NULL, maxp=min(nrow(x),ncol(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitobj}{
Result of a call to \code{lars} or  \code{lars.en} or  \code{lars.glm}
}
  \item{x}{
N by p matrix of predictors
     }
  \item{y}{
N-vector of outcome values
}
  \item{sigma.est}{
Estimate of error standard deviation. If a numerical value, that value if used. If "full"
the  (square root)  of the mean squared residual from the full model is used.
}
  \item{status}{
Optional N-vector of censoring indicators for Cox Proportional hazards model. 1=failed; 0=censored.
}
\item{maxp}{
Optional limit for number of steps to be analyzed.
}}
\details{
This function  computes covariance test for inference in  adaptive linear modelling, for lasso (least angle
regression) in the Gaussian case, binomial/logistic and Cox proportional hazards survival models. It estimates p-values for each predictor entered, that
account for the adpative nature of the fitting.
}
\value{
  \item{results }{Table of covariance test values and p-values, for each predictor entered}
  \item{sigma}{Estimate of sigma used}
  \item{null.dist}{Null distribution used to compute p-values}
}
\references{
A significance test for the lasso (2013). Lockhart, R.,  Taylor, J.,  Tibshirani (Ryan) and  Tibshirani (Robert)
}
\author{
Rob Tibshirani
}
\note{
This function requires the \code{lars} R library (for the Gaussian case), and the \code{glmpath}
function for the logistic and Cox model.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
lars, lars.en, lars.glm
}
\examples{


set.seed(1234)
x=matrix(rnorm(100*10),ncol=10)
mx=colMeans(x)
sx=sqrt(apply(x,2,var))
x=scale(x,mx,sx)
beta=c(4,rep(0,9))
y=x\%*\%beta+.4*rnorm(100)

#Gaussian
a=lars(x,y)
covTest(a,x,y)

#EN
a=lars.en(x,y,lambda2=1)
covTest(a,x,y)

#logistic
y=1*(y>0)

a=lars.glm(x,y,family="binomial")
covTest(a,x,y)


# Cox model
#y=6*x[,2]+rnorm(100)+10
#status=sample(c(0,1),size=length(y),replace=TRUE)
#a=lars.glm(x,y,status=status,family="cox")
#covTest(a,x,y,status=status)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{lasso}
\keyword{least angle regression} 
\keyword{adaptive fitting} 

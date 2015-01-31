\name{covTest-package}
\alias{covTest-package}
\alias{covTest-packarghe}
\docType{package}
\title{
Compute the covariance test for adaptive linear modelling
}
\description{
Compute the covariance test significance testing in  adaptive linear modelling.
Can be used with LARS (lasso) for linear models, elastic net, binomial
and Cox survival model.
This package should be considered EXPERIMENTAL. The background paper is not yet published
and rigorous theory does not yet exist for the logistic and Cox models.

}
\details{
\tabular{ll}{
Package: \tab covTest\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-01-08\cr
License: \tab GPL-2\cr
}
Very simple to use. Takes output from one of \code{lars}, \code{lars.en},\code{lars.glm} and compute covariance test and p-values.
Requires \code{lars} and \code{glmpath} packages. \code{lars.en} and \code{lars.glm}
are included in this package.
Functions are:
\code{covTest}
\code{lars.en}
\code{lars.glm}
\code{predict.lars.en}
\code{predict.lars.glm}
}
\author{Rob Tibshirani tibs@stanford.edu
}
\references{
A significance test for the lasso (2013). Lockhart, R.,  Taylor, J.,  Tibshirani (Ryan) and  Tibshirani (Robert)
}
\keyword{ LARS }
\seealso{covTest, lars.glm, lars.en
}
\examples{

x=matrix(rnorm(100*10),ncol=10)
x=scale(x,TRUE,TRUE)/sqrt(99)

#Gaussian
beta=c(4,rep(0,9))
y=x\%*\%beta+.4*rnorm(100)

a=lars(x,y)
covTest(a,x,y)

#Elastic net
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


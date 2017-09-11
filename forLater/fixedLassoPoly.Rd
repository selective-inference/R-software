\name{fixedLassoPoly}
\alias{fixedLassoPoly}

\title{
Compute polyhedral constraints for a LASSO problem with
a fixed value of lambda.
}
\description{
Compute polyhedral representation of the selection region of Lee et al. (2016).
By construction, y should satisfy A %*% y elementwise less then or equal b.
}
\usage{
fixedLassoPoly(X, y, lambda, beta, active, inactive = FALSE)
}
\arguments{
\item{X}{
Design matrix of LASSO problem.
}
\item{y}{
Response of LASSO problem.
}
\item{lambda}{
Value of regularization parameter.
}      
\item{beta}{
Solution of LASSO problem with regularization parameter set to lambda.
}      
\item{active}{
Active set of the LASSO problem as a boolean vector. Should correspond
to the non-zeros of beta.
} 
\item{inactive}{
Form the inactive constraints as well?
} 
}
\details{
This function computes
the  polyhedral representation of the selection region of Lee et al. (2016).
}

\value{  
\item{A}{Linear part of the affine inequalities.}
\item{b}{RHS offset the affine inequalities.}
}

\references{
Jason Lee, Dennis Sun, Yuekai Sun, and Jonathan Taylor (2016). 
Exact post-selection inference, with application to the lasso. Annals of Statistics, 44(3), 907-927.

Jonathan Taylor and Robert Tibshirani (2017) Post-selection inference for math L1-penalized likelihood models.
Canadian Journal of Statistics, xx, 1-21. (Volume still not posted)
}
\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid}

\examples{

set.seed(43)
n = 50
p = 10
sigma = 1

x = matrix(rnorm(n*p),n,p)
x = scale(x,TRUE,TRUE)

beta = c(3,2,rep(0,p-2))
y = x\%*\%beta + sigma*rnorm(n)

# first run glmnet
gfit = glmnet(x,y,standardize=FALSE)

# extract coef for a given lambda; note the 1/n factor!
# (and we don't save the intercept term)
lambda = .8
beta = coef(gfit, s=lambda/n, exact=TRUE)[-1]
active = (beta != 0)

fixedLassoPoly(x, y, lambda, beta, active)
fixedLassoPoly(x, y, lambda, beta, active, inactive=TRUE)

}
 
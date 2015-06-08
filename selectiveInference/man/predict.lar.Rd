\name{predict.lar}
\alias{predict.lar}
\title{
Prediction from a LAR fit
}
\description{
Make predictions or extract coefficients from a fitted lar model
}
\usage{
\method{predict}{lar}(object, newx, s, type = c("fit", "coefficients"), mode = c("step", "fraction", "norm", "lambda"), ...)
}

\arguments{
  \item{object}{
Object returned by a call to \code{lar} function (Not the \code{lars} function)
}
  \item{newx}{
If type="fit", then newx should be the x values at which the
          fit is required. If type="coefficients", then newx can be
          omitted.
}
  \item{s}{
        a value, or vector of values, indexing the path. Its values
          depends on the mode= argument. By default (mode="step"), s
          should take on values between 0 and p (e.g., a step of 1.3
          means .3 of the way between step 1 and 2.)
}
  \item{type}{
If type="fit", predict returns the fitted values. If
          type="coefficients", predict returns the coefficients.
          Abbreviations allowed.
}
  \item{mode}{
Mode="step" means the s= argument indexes the lars step
          number, and the coefficients will be returned corresponding
          to the values corresponding to step s. If mode="fraction",
          then s should be a number between 0 and 1, and it refers to
          the ratio of the L1 norm of the coefficient vector, relative
          to the norm at the full LS solution. Mode="norm" means s
          refers to the L1 norm of the coefficient vector.
          Mode="lambda" uses the lasso regularization parameter for s;
          for other models it is the maximal correlation (does not make
          sense for lars/stepwise models). Abbreviations allowed.
}
  \item{\dots}{
Any arguments for predict.lar should work for coef.lars.
        }
}
\details{
LARS is described in detail in Efron, Hastie, Johnstone and
     Tibshirani (2002).
}
\value{
Either a vector/matrix of fitted values, or a vector/matrix of
     coefficients.
}
\references{
Efron, Hastie, Johnstone and Tibshirani (2002) "Least Angle
     Regression" (with discussion) _Annals of Statistics_; see also
  
     http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf

     Hastie, Tibshirani and Friedman (2002) Elements of Statistical
     Learning, Springer, NY.

      Jonathan Taylor, Richard Lockhart,  Ryan Tibshirani and Rob Tibshirani (2015) 
    Exact Post-selection Inference for Forward Stepwise and Least Angle Regression.
arXiv:1401.3889
}
\author{
Joshua Loftus, Stephen Reid, Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
}

\examples{
#NOT RUN
#set.seed(33)
#n=200
#p=20
#sigma=1
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,T)/sqrt(n-1)
#generate y
#beta=c(3,-2,rep(0,p-2))
#beta=c(rep(3,10),rep(0,p-10))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#larfit=lar(x,y)
#fit=predict.lar(larfit,x,type="fit")
}

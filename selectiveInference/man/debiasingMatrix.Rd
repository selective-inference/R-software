\name{debiasingMatrix}
\alias{debiasingMatrix}
\title{
Find an approximate inverse of a non-negative definite matrix.
}
\description{
Find some rows of an approximate inverse of a non-negative definite 
symmetric matrix by solving optimization problem described 
in Javanmard and Montanari (2013). Can be used for approximate 
Newton step from some consistent estimator (such as the LASSO)
to find a debiased solution.
}
\usage{
debiasingMatrix(Xinfo, 
                is_wide,			
                nsample, 
                rows, 
		verbose=FALSE, 
		bound=NULL,            
   		linesearch=TRUE,    
   		scaling_factor=1.5, 
		max_active=NULL,    
		max_try=10,         
		warn_kkt=FALSE,     
		max_iter=50,       
		kkt_stop=TRUE,
		parameter_stop=TRUE,
		objective_stop=TRUE,
                kkt_tol=1.e-4,      
		parameter_tol=1.e-4,
		objective_tol=1.e-4)
}
\arguments{
\item{Xinfo}{
Either a non-negative definite matrix S=t(X) %*% X / n or X itself. If 
is_wide is TRUE, then Xinfo should be X, otherwise it should be S.
}
\item{is_wide}{
Are we solving for rows of the debiasing matrix assuming it is 
a wide matrix so that Xinfo=X and the non-negative definite
matrix of interest is t(X) %*% X / nrow(X).
}
\item{nsample}{
Number of samples used in forming the cross-covariance matrix.
Used for default value of the bound parameter.
}
\item{rows}{
Which rows of the approximate inverse to compute.
}      
\item{verbose}{
Print out progress as rows are being computed.
} 
\item{bound}{
Initial bound parameter for each row. Will be changed
if linesearch is TRUE.
} 
\item{linesearch}{
Run a line search to find as small as possible a bound parameter for each row?
}
\item{scaling_factor}{
In the linesearch, the bound parameter is either multiplied or divided by this
factor at each step.
}
\item{max_active}{
How large an active set to consider in solving the problem with coordinate descent.
Defaults to max(50, 0.3*nsample).
}
\item{max_try}{
How many tries in the linesearch.
}
\item{warn_kkt}{
Warn if the problem does not seem to be feasible after running the coordinate
descent algorithm.
}
\item{max_iter}{
How many full iterations to run of the coordinate descent for each
value of the bound parameter.
}
\item{kkt_stop}{
If TRUE, check to stop coordinate descent when KKT conditions are approximately satisfied.
}
\item{parameter_stop}{
If TRUE, check to stop coordinate descent based on relative convergence of parameter vector,
checked at geometrically spaced iterations 2^k.
}
\item{objective_stop}{
If TRUE, check to stop coordinate descent based on relative decrease of objective value,
checked at geometrically spaced iterations 2^k.
}
\item{kkt_tol}{
Tolerance value for assessing whether KKT conditions for solving the
dual problem and feasibility of the original problem.
}
\item{parameter_tol}{
Tolerance value for assessing convergence of the problem using relative
convergence of the parameter.
}
\item{objective_tol}{
Tolerance value for assessing convergence of the problem using relative
decrease of the objective.
}
}
\details{
This function computes an approximate inverse
as described in Javanmard and Montanari (2013), specifically
display (4). The problem is solved by considering a dual
problem which has an objective similar to a LASSO problem and is solvable
by coordinate descent. For some values of bound the original
problem may not be feasible, in which case the dual problem has no solution. 
An attempt to detect this is made by stopping when the active set grows quite
large, determined by max_active.
}

\value{  
\item{M}{Rows of approximate inverse of Sigma.}
}

\references{
Adel Javanmard and Andrea Montanari (2013). 
Confidence Intervals and Hypothesis Testing for High-Dimensional Regression. Arxiv: 1306.3171
}
\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid}

\examples{

set.seed(10)
n = 50
p = 100
X = matrix(rnorm(n * p), n, p)
S = t(X) \%*\% X / n
M = debiasingMatrix(S, FALSE, n, c(1,3,5))
M2 = debiasingMatrix(X, TRUE, n, c(1,3,5))
max(M - M2)
}
 
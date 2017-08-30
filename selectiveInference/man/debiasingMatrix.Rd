\name{debiasingMatrix}
\alias{debiasingMatrix}
\title{
Find an approximate inverse of a non-negative definite matrix.
}
\description{
Find some rows of an approximate inverse of a non-negative definite 
symmetric matrix by solving optimization problem described 
in Javanmard and Montanari (2013).
}
\usage{
debiasingMatrix(Sigma, 
                nsample, 
                rows, 
		verbose=FALSE, 
		mu=NULL,            
   		linesearch=TRUE,    
   		scaling_factor=1.5, 
		max_active=NULL,    
		max_try=10,         
		warn_kkt=FALSE,     
		max_iter=100,       
                kkt_tol=1.e-4,      
		objective_tol=1.e-8)
}
\arguments{
\item{Sigma}{
A symmetric non-negative definite matrix, often a cross-covariance matrix.
}
\item{nsample}{
Number of samples used in forming the cross-covariance matrix.
Used for default value of the bound parameter mu.
}
\item{rows}{
Which rows of the approximate inverse to compute.
}      
\item{verbose}{
Print out progress as rows are being computed.
} 
\item{mu}{
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
\item{kkt_tol}{
Tolerance value for assessing whether KKT conditions for solving the
dual problem and feasibility of the original problem.
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
by coordinate descent. For some values of mu the original
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
M = debiasingMatrix(S, n, c(1,3,5))

}
 
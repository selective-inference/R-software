\name{randomizedLasso}
\alias{randomizedLasso}

\title{
Inference for the randomized lasso, with a fixed lambda 
}
\description{
Solve a randomly perturbed LASSO problem.
}
\usage{
randomizedLasso(X, 
                y, 
                lam, 
                family=c("gaussian", "binomial"),
                noise_scale=NULL, 
                ridge_term=NULL, 
                max_iter=100,       
                kkt_tol=1.e-4,      
                parameter_tol=1.e-8,
                objective_tol=1.e-8,
                objective_stop=FALSE,
                kkt_stop=TRUE,
                parameter_stop=TRUE)
}
\arguments{
  \item{X}{
Matrix of predictors (n by p); 
}
  \item{y}{
Vector of outcomes (length n)
}
  \item{lam}{
Value of lambda used to compute beta. See the above warning
 Be careful! This function uses the "standard" lasso objective
  \deqn{
    1/2 \|y - x \beta\|_2^2 + \lambda \|\beta\|_1.
  }
 In contrast, glmnet multiplies the first term by a factor of 1/n.
 So after running glmnet, to extract the beta corresponding to a value lambda, 
 you need to use \code{beta = coef(obj, s=lambda/n)[-1]},
 where obj is the object returned by glmnet (and [-1] removes the intercept,
 which glmnet always puts in the first component)
} 
\item{family}{
Response type: "gaussian" (default), "binomial".
}
\item{noise_scale}{
Scale of Gaussian noise added to objective. Default is 
0.5 * sd(y) times the sqrt of the mean of the trace of X^TX.
}
\item{ridge_term}{
A small "elastic net" or ridge penalty is added to ensure
the randomized problem has a solution. 
0.5 * sd(y) times the sqrt of the mean of the trace of X^TX divided by
sqrt(n).
}
\item{max_iter}{
How many rounds of updates used of coordinate descent in solving randomized
LASSO.
}
\item{kkt_tol}{
Tolerance for checking convergence based on KKT conditions.
}
\item{parameter_tol}{
Tolerance for checking convergence based on convergence
of parameters.
}
\item{objective_tol}{
Tolerance for checking convergence based on convergence
of objective value.
}
\item{kkt_stop}{
Should we use KKT check to determine when to stop?
}
\item{parameter_stop}{
Should we use convergence of parameters to determine when to stop?
}
\item{objective_stop}{
Should we use convergence of objective value to determine when to stop?
}
}

\details{For \code{family="gaussian"} this function uses the "standard" lasso objective
  \deqn{
    1/2 \|y - x \beta\|_2^2 + \lambda \|\beta\|_1
  }
and adds a term 
  \deqn{
    - \omega^T\beta + \frac{\epsilon}{2} \|\beta\|^2_2
  }
where omega is drawn from IID normals with standard deviation
\code{noise_scale} and epsilon given by \code{ridge_term}. 
See below for default values of \code{noise_scale} and \code{ridge_term}.

For \code{family="binomial"}, the squared error loss is replaced by the
negative of the logistic log-likelihood.
}
\value{  
\item{X}{Design matrix.}
\item{y}{Response vector.}
\item{lam}{Vector of penalty parameters.}
\item{family}{Family: "gaussian" or "binomial".}
\item{active_set}{Set of non-zero coefficients in randomized solution that were penalized. Integers from 1:p.}
\item{inactive_set}{Set of zero coefficients in randomized solution. Integers from 1:p.}
\item{unpenalized_set}{Set of non-zero coefficients in randomized solution that were not penalized. Integers from 1:p.}
\item{sign_soln}{The sign pattern of the randomized solution.}
\item{full_law}{List describing sampling parameters for conditional law of all optimization variables given the data in the LASSO problem.}
\item{conditional_law}{List describing sampling parameters for conditional law of only the scaling variables given the data and the observed subgradient in the LASSO problem.}
\item{internal_transform}{Affine transformation describing relationship between internal representation of the data and the data compontent of score of the likelihood at the unregularized MLE based on the sign_vector (a.k.a. relaxed LASSO).}
\item{observed_raw}{Data component of the score at the unregularized MLE.}
\item{noise_scale}{SD of Gaussian noise used to draw the perturbed objective.}
\item{soln}{The randomized solution. Inference is made conditional on its sign vector (so no more snooping of this value is formally permitted.) 
If \code{condition_subgrad == TRUE} when sampling, then we may snoop on the observed subgradient.}
\item{perturb}{The random vector in the linear term added to the objective.}
}

\references{
Xiaoying Tian, and Jonathan Taylor (2015).
Selective inference with a randomized response. arxiv.org:1507.06739

Xiaoying Tian, Snigdha Panigrahi, Jelena Markovic, Nan Bi and Jonathan Taylor (2016).
Selective inference after solving a convex problem. 
arxiv:1609.05609

}
\author{Jelena Markovic, Jonathan Taylor}

\examples{
set.seed(43)
n = 50
p = 10
sigma = 0.2
lam = 0.5

X = matrix(rnorm(n*p), n, p)
X = scale(X, TRUE, TRUE) / sqrt(n-1)

beta = c(3,2,rep(0,p-2))
y = X\%*\%beta + sigma*rnorm(n)

result = randomizedLasso(X, y, lam)

}
 

\name{fixedLassoInf}
\alias{fixedLassoInf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inference for fixed lambda lasso
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fixedLassoInf(x, y, bhat, lam, sigma, alpha = 0.1, trace = F, compute.ci = F, tol = 1e-05, one.sided = TRUE, nsigma = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{bhat}{
%%     ~~Describe \code{bhat} here~~
}
  \item{lam}{
%%     ~~Describe \code{lam} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
}
  \item{compute.ci}{
%%     ~~Describe \code{compute.ci} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{one.sided}{
%%     ~~Describe \code{one.sided} here~~
}
  \item{nsigma}{
%%     ~~Describe \code{nsigma} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, bhat, lam, sigma, alpha = 0.1, trace = F, compute.ci = F, 
    tol = 1e-05, one.sided = TRUE, nsigma = 10) 
{
    this.call = match.call()
    junk = tf.jonab(y, x, bhat, lam, tol = tol)
    n = length(y)
    a = junk$A
    b = junk$b
    e = which(bhat != 0)
    xe = x[, e, drop = F]
    pp = length(e)
    pv = vmall = vpall = rep(NA, pp)
    ci = miscov = matrix(NA, pp, 2)
    etaall = matrix(NA, nrow = pp, ncol = n)
    SMALL = 1e-07
    for (k in 1:pp) {
        if (trace) 
            cat(k, fill = T)
        eta = ginv(t(xe))[, k]
        bhat0 = sum(eta * y)
        if (one.sided) 
            eta = eta * sign(bhat0)
        etaall[k, ] = eta
        vs = tf.jonvs(y, a, b, eta)
        vpp = vs$vp
        vmm = vs$vm
        vmall[k] = vmm
        vpall[k] = vpp
        tt = sum(eta * y)
        sigma.eta = sigma * sqrt(sum(eta^2))
        u = 0
        pv[k] = 1 - ptruncnorm(tt, vmm, vpp, u, sigma.eta)
        if (!one.sided) 
            pv[k] = 2 * min(pv[k], 1 - pv[k])
        if (compute.ci) {
            junk = selection.int(y, eta, sigma^2, vs, alpha, 
                nsigma = nsigma)
            ci[k, ] = junk$ci
            miscov[k, ] = junk$miscov
        }
    }
    out = list(lam = lam, eta = etaall, vm = vmall, vp = vpall, 
        pv = pv, ci = ci, miscov = miscov, pred = e, alpha = alpha, 
        lambda = lam)
    class(out) = "fixedLassoInf"
    out$call = this.call
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

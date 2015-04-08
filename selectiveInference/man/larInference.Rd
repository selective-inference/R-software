\name{larInference}
\alias{larInference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selective inference for least angle regression (LAR)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
larInference(x, y, larfit, sigma, compute.ci = TRUE, alpha = 0.1, one.sided = TRUE, nsigma = 10, maxp = min(nrow(x), ncol(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{larfit}{
%%     ~~Describe \code{larfit} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{compute.ci}{
%%     ~~Describe \code{compute.ci} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{one.sided}{
%%     ~~Describe \code{one.sided} here~~
}
  \item{nsigma}{
%%     ~~Describe \code{nsigma} here~~
}
  \item{maxp}{
%%     ~~Describe \code{maxp} here~~
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
function (x, y, larfit, sigma, compute.ci = TRUE, alpha = 0.1, 
    one.sided = TRUE, nsigma = 10, maxp = min(nrow(x), ncol(x))) 
{
    this.call = match.call()
    SMALL = 1e-06
    p = ncol(x)
    nk = larfit$nk
    vmm = vpp = pv = sigma.eta = rep(NA, maxp)
    ci = miscov = matrix(NA, maxp, 2)
    for (k in 1:maxp) {
        mod = larfit$act[1:k]
        temp = (solve(t(x[, mod, drop = F]) \%*\% x[, mod, drop = F]) \%*\% 
            t(x[, mod, drop = F]))
        temp = temp[nrow(temp), ]
        bhat = sum(temp * y)
        eta = as.vector(temp)
        if (one.sided) 
            eta = eta * sign(bhat)
        tt = sum(eta * y)
        A = -larfit$Gam[1:nk[k], ]
        pp = nrow(A)
        b = rep(0, pp)
        alp = as.vector(A \%*\% eta/sum(eta^2))
        alp[abs(alp) < SMALL] = 0
        vp = rep(Inf, pp)
        vm = rep(-Inf, pp)
        for (jj in 1:pp) {
            if (alp[jj] < 0) 
                vm[jj] = (b[jj] - (A \%*\% y)[jj] + alp[jj] * sum(eta * 
                  y))/alp[jj]
            if (alp[jj] > 0) 
                vp[jj] = (b[jj] - (A \%*\% y)[jj] + alp[jj] * sum(eta * 
                  y))/alp[jj]
        }
        vmm[k] = max(vm, na.rm = T)
        vpp[k] = min(vp, na.rm = T)
        sigma.eta[k] = sigma * sqrt(sum(eta^2))
        u = 0
        pv[k] = 1 - ptruncnorm(tt, vmm[k], vpp[k], u, sigma.eta[k])
        if (!one.sided) 
            pv[k] = 2 * min(pv[k], 1 - pv[k])
        if (compute.ci) {
            vs = list(vm = vmm[k], vp = vpp[k])
            junk = selection.int(y, eta, sigma, vs, alpha, nsigma = nsigma)
            ci[k, ] = junk$ci
            miscov[k, ] = junk$miscov
        }
    }
    pv.spacing = spacing.pval.asymp.list(y, larfit, maxp, sigma = sigma)
    junk = covtest(larfit, x, y, sigma, maxp)
    pv.cov = 1 - pexp(junk, 1)
    out = list(pv = pv, ci = ci, miscov = miscov, vmm = vmm, 
        vpp = vpp, sigma.eta = sigma.eta, alpha = alpha, act = larfit$act, 
        pv.spacing = pv.spacing, pv.cov = pv.cov)
    out$call = this.call
    class(out) = "larInference"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

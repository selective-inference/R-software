\name{selectiveInference-internal}
\title{Internal PMA functions}
 \alias{selection.int}
\alias{compute.vmvp}
 \alias{myptruncnorm}
 \alias{plot.lar} 
\alias{print.fixedLassoInf}
  \alias{print.forwardStep} 
\alias{print.forwardStepInf}
\alias{print.lar} 
\alias{print.larInf}
  \alias{print.manyMeans}


\description{Internal selectiveInference functions}
\usage{
selection.int(y, eta, sigma, vs, alpha, del = 1e-04, flip = F)
compute.vmvp(y, eta, A, b, pp, SMALL = 1e-07)
myptruncnorm(z, a, b, mean = 0, sd = 1)
\method{plot}{lar}(x, xvar = c("norm", "step", "lambda"), breaks = TRUE,
                 omit.zeros = TRUE, ...)
\method{plot}{larInf}{x, ...}
\method{print}{fixedLassoInf}(x, digits = max(3, getOption("digits") - 3), ...)
\method{print}{forwardStep}(x, digits = max(3, getOption("digits") - 3), ...)
\method{print}{forwardStepInf}(x, digits = max(3, getOption("digits") - 3), ...)
\method{print}{lar}(x,...)
\method{print}{larInf}(x, digits = max(3, getOption("digits") - 3), ...)
\method{print}{manyMeans}(x,...)
}
\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Max G'Sell, Stephen Reid}

\keyword{internal}



check_inequalities <- function(fit, x, y, index, k, R, tol = 1e-15) {

  U <- R/sqrt(sum(R^2))
  L <- lapply(1:length(fit$projections), function(l) {

      Q <- fit$projections[[l]]
      # The quadratic form corresponding to
      # (t*U + Z)^T %*% Q %*% (t*U + Z) \geq 0
      # we find the roots in t, if there are any
      # and return the interval of potential t
      kl <- ifelse(k > 0, k * sum(diag(Q)), 0)
      if (t(y) %*% Q %*% y - kl < -.Machine$double.eps) {
        print(paste("Problematic projection:", l))
        stop("Observation does not belong to selection region")
      }

      Z <- y - R
      a = t(U) %*% Q %*% U
      b = 2 * t(U) %*% Q %*% Z
      c = t(Z) %*% Q %*% Z - kl
      disc <- b^2 - 4*a*c
      b2a <- -b/(2*a)

      if (disc > tol) {
        # Real roots
        pm <- sqrt(disc)/(2*a)
        endpoints <- c(b2a - pm, b2a + pm)
        
      } else {
            
        # No real roots

        #if (a > 0) {
          # Quadratic form always positive
          return(Intervals(c(0,Inf)))
        #} else { #### Assume this never happens ####
          # Quadratic form always negative
        #  stop("Infeasible!")
        #}
      }
          
      if (a > tol) {
        # Parabola opens upward
            
        if (min(endpoints) > 0) {
          # Both roots positive, union of intervals
          return(Intervals(rbind(c(0, min(endpoints)), c(max(endpoints), Inf))))
            
        } else {
          # At least one negative root
          return(Intervals(c(max(0, max(endpoints)), Inf)))
        }
         
      } else {
        if (a < -tol) {
          # Parabola opens downward
          if (max(endpoints) < 0) {
             
            # Positive quadratic form only when t negative
            stop("Error: infeasible")
          } else {
              
            # Part which is positive
            return(Intervals(c(max(0, min(endpoints)), max(endpoints))))
          }
        } else {
          # a is too close to 0, quadratic is actually linear
          if (abs(b) > tol) {
            if (b > 0) {
              return(Intervals(c(max(0, -c/b), Inf)))
            } else {
              if (-c/b < 0) stop("Error: infeasible linear equation")
              return(Intervals(c(0, -c/b)))
            }
          } else {
            warning("Ill-conditioned quadratic")
            return(Intervals(c(0, Inf)))
          }
            
        }
      }
    })
    # L is now a list of lists of intervals
    return(L)
}

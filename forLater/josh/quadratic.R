
check_inequalities <- function(obj, x, y, index, k, TC, R, Ugtilde, tol = 1e-15) {

  eta <- Ugtilde %*% R / TC

  L <- lapply(1:length(obj$action), function(step) {
      
    Ug <- obj$maxprojs[[step]]
  
    lapply(1:length(obj$projections[[step]]), function(l) {

      Uh <- obj$projections[[step]][[l]]
      # The quadratic form corresponding to
      # (t*U + Z)^T %*% Q %*% (t*U + Z) \geq 0
      # we find the roots in t, if there are any
      # and return the interval of potential t
      kl <- ifelse(k > 0, k * ncol(Uh), 0)
      Uhy <- t(Uh) %*% y
      if (sum(Uhy^2) - kl < -.Machine$double.eps) {
        print(paste("Problematic projection:", l))
        stop("Observation does not belong to selection region")
      }

      Z <- y - eta * TC
      Uheta <- t(Uh) %*% eta
      Ugeta <- t(Ug) %*% eta
      UhZ <- t(Uh) %*% Z
      UgZ <- t(Ug) %*% Z
      a = sum(Ugeta^2) - sum(Uheta^2)
      b = 2 * (t(Ugeta) %*% UgZ - t(Uheta) %*% UhZ)
      c = sum(UgZ^2) - sum(UhZ^2) - kl 

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
          return(Intervals(c(-Inf,0)))
        #} else { #### Assume this never happens ####
          # Quadratic form always negative
        #  stop("Infeasible!")
        #}
      }
          
      if (a > tol) {
        # Parabola opens upward
        if (min(endpoints) > 0) {
          # Both roots positive, union of intervals
          return(Intervals(rbind(c(-Inf,0), c(min(endpoints), max(endpoints)))))
        } else {
          # At least one negative root
          return(Intervals(c(-Inf, max(0, max(endpoints)))))
        }
      } else {
        if (a < -tol) {
          # Parabola opens downward
          if (max(endpoints) < 0) {
             
            # Positive quadratic form only when t negative
            stop("Error: infeasible")
          } else {
            # Part which is positive
            return(Intervals(rbind(c(-Inf, max(0, min(endpoints))), c(max(endpoints), Inf))))
          }
        } else {
          # a is too close to 0, quadratic is actually linear
          if (abs(b) > tol) {
            if (b > 0) {
              return(Intervals(c(-Inf, max(0, -c/b))))
            } else {
              if (-c/b < 0) stop("Error: infeasible linear equation")
              return(Intervals(rbind(c(-Inf, 0), c(-c/b, Inf))))
            }
          } else {
            warning("Ill-conditioned quadratic")
            return(Intervals(c(-Inf,0)))
          }
        }
      }
    })
    # LL is a list of intervals
  })
  # L is now a list of lists of intervals
  return(unlist(L, recursive = FALSE, use.names = FALSE))
}


roots_to_checkpoints <- function(roots) {
    checkpoints <- unique(sort(c(0, roots)))
    return(c(0, (checkpoints + c(checkpoints[-1], 2 + checkpoints[length(checkpoints)]))/2))
}

roots_to_partition <- function(roots) {
    checkpoints <- unique(sort(c(0, roots)))
    return(list(endpoints = c(checkpoints, Inf), midpoints = (checkpoints + c(checkpoints[-1], 2 + checkpoints[length(checkpoints)]))/2))
}

# Tchi_roots <- function(Q, a, b, )

TF_roots <- function(Q, a, b, Vdelta, V2, z, C, r, tol = 1e-14) {

    # z = y - R1
    VdeltaQ <- t(Vdelta) %*% Q
    V2Q <- t(V2) %*% Q
    x11 <- VdeltaQ %*% Vdelta
    x12 <- 2 * VdeltaQ %*% V2
    x22 <- V2Q %*% V2
    x1 <- 2 * VdeltaQ %*% z + t(a) %*% Vdelta
    x2 <- 2 * V2Q %*% z + t(a) %*% V2
    x0 <- t(z) %*% Q %*% z + t(a) %*% z + b

    g1 <- function(t) r*sqrt(C*t/(1+C*t))
    g2 <- function(t) r/sqrt(1+C*t)
    I <- function(t) x11*g1(t)^2 + x12*g1(t)*g2(t) + x22*g2(t)^2 + x1*g1(t) + x2*g2(t) + x0

    z4 <- r*complex(real = -x11 + x22, imaginary = -x12)/2
    z3 <- complex(real = x2, imaginary = -x1)
    z2 <- complex(real = r*x11+r*x22+2*x0/r)
    z1 <- Conj(z3)
    z0 <- Conj(z4)
    zcoefs <- r*c(z0, z1, z2, z3, z4)/2
    croots <- polyroot(zcoefs)
    thetas <- Arg(croots)
    modinds <- Mod(croots) <= 1 + tol & Mod(croots) >= 1 - tol
    angleinds <- thetas >=0 & thetas <= pi/2
    roots <- unique(thetas[modinds * angleinds])
    troots <- tan(roots)^2/C
    
    if (length(roots) == 0) {
        return(list(intervals = Intervals(c(0,Inf)), I=I))
    }
    
    checkpoints <- roots_to_checkpoints(troots)
    signs <- sign(I(checkpoints))
    diffs <- c(0, diff(signs))
    changeinds <- which(diffs != 0)
    
    if (length(changeinds) > 0) {
        
        roots <- unlist(lapply(changeinds, function(ind) {
            uniroot(I, lower = checkpoints[ind-1], upper = checkpoints[ind])$root
        }))
        partition <- roots_to_partition(roots)
        positive <- which(I(partition$midpoints) > 0)
        
        intervals <- matrix(NA, ncol=2)
        for (i in 1:length(positive)) {
            ind <- positive[i]
            if ((i > 1) && (ind == positive[i-1] + 1)) {
                intervals[nrow(intervals), 2] <- partition$endpoints[ind+1]
            } else {
                intervals <- rbind(intervals, c(partition$endpoints[ind], partition$endpoints[ind+1]))
            }
        }

        return(list(intervals = Intervals(intervals[-1,]), I=I))
    }
    
    return(list(intervals = Intervals(c(0,Inf)), I=I))
}


## # test
## for(i in 1:10) {
## n <- 100
## p <- 90
## y <- rnorm(n)
## x <- rnorm(n)
## X2 <- matrix(rnorm((p-1)*n),ncol=p-1)
## X <- cbind(x,X2)
## Psub <- x %*% ginv(x)
## Pfull <- X %*% ginv(X)
## Pz <- diag(rep(1,n)) - Psub
## PM <- diag(rep(1,n)) - Pfull
## R1 <- Pz %*% y
## R2 <- PM %*% y
## z <- y - R1
## C <- sum(diag(Pfull-Psub))/sum(diag(Psub))
## norm2R1 <- sum(R1^2)
## norm2R2 <- sum(R2^2)
## TF <- (norm2R1-norm2R2)/(C*norm2R2)
## r = sqrt(norm2R2)
## Vdelta <- (R1-R2)/sqrt(norm2R1-norm2R2)
## V2 <- R2/r
## Q <- diag(runif(n))
## a <- rnorm(n)
## b <- rnorm(1)
## froots <- F_roots(Q,a,b,Vdelta,V2,z,C,r)
## print(froots$intervals)
## }
## #        t <- seq(from = 0, to = 50, length.out = 10000)
## #        plot(t,froots$I(t), type = "l")
## #        points(roots,I(roots),col="red")
## #        abline(h=0)

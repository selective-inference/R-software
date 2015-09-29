
interval_groupfs <- function(obj, TC, R, eta, Ugtilde, tol = 1e-15) {

  Z <- obj$y - eta * TC
  n <- nrow(obj$x)

  L <- lapply(1:length(obj$action), function(s) {

    Ug <- obj$maxprojs[[s]]
    dfg <- ncol(Ug)

    if (s > 1) {
        etas <- obj$cumprojs[[s-1]] %*% eta
        Zs <- obj$cumprojs[[s-1]] %*% Z
    } else {
        etas <- eta
        Zs <- Z
    }

    num.projs <- length(obj$projections[[s]])
    if (num.projs == 0) {
        return(list(Intervals(c(-Inf,0))))
    } else {
      lapply(1:num.projs, function(l) {

      Uh <- obj$projections[[s]][[l]]
      dfh <- ncol(Uh)
      # The quadratic form corresponding to
      # (t*U + Z)^T %*% Q %*% (t*U + Z) \geq 0
      # we find the roots in t, if there are any
      # and return the interval of potential t

      Uheta <- t(Uh) %*% etas
      Ugeta <- t(Ug) %*% etas
      UhZ <- t(Uh) %*% Zs
      UgZ <- t(Ug) %*% Zs
      etasZs <- t(etas) %*% Zs
      peng <- obj$maxpens[[s]]
      penh <- obj$aicpens[[s]][[l]]
      pendiff <- peng-penh
      if (is.null(obj$sigma)) {
          A <- sum(Ugeta^2) * peng - sum(Uheta^2) * penh - sum(etas^2) * pendiff
          B <- 2 * as.numeric(t(Ugeta) %*% UgZ * peng - t(Uheta) %*% UhZ * penh - etasZs * pendiff)
          C <- sum(UgZ^2) * peng - sum(UhZ^2) * penh - sum(Zs^2) * pendiff
      } else {
          # Check this
          A <- sum(Ugeta^2) - sum(Uheta^2)
          B <- 2 * as.numeric(t(Ugeta) %*% UgZ - t(Uheta) %*% UhZ)
          C <- sum(UgZ^2) - sum(UhZ^2) - pendiff
      }

      disc <- B^2 - 4*A*C
      b2a <- -B/(2*A)

      if (disc > tol) {
        # Real roots
        pm <- sqrt(disc)/(2*A)
        endpoints <- sort(c(b2a - pm, b2a + pm))

      } else {

        # No real roots
        if (A > -tol) {
          # Quadratic form always positive
          return(Intervals(c(-Inf,0)))
        } else {
          # Quadratic form always negative
          stop(paste("Empty TC support is infeasible", s, "-", l))
        }
      }

      if (A > tol) {
        # Parabola opens upward
        if (min(endpoints) > 0) {
          # Both roots positive, union of intervals
          return(Intervals(rbind(c(-Inf,0), endpoints)))
        } else {
          # At least one negative root
          return(Intervals(c(-Inf, max(0, endpoints[2]))))
        }
      } else {
        if (A < -tol) {
          # Parabola opens downward
          if (endpoints[2] < 0) {
            # Positive quadratic form only when t negative
            stop(paste("Negative TC support is infeasible", s, "-", l))
          } else {
            # Part which is positive
            if (endpoints[1] > 0) {
                return(Intervals(rbind(c(-Inf, endpoints[1]), c(endpoints[2], Inf))))
            } else {
                return(Intervals(c(endpoints[2], Inf)))
            }
          }
        } else {
          # a is too close to 0, quadratic is actually linear
          if (abs(B) > tol) {
            if (B > 0) {
              return(Intervals(c(-Inf, max(0, -C/B))))
            } else {
              if (-C/B < 0) stop("Error: infeasible linear equation")
              return(Intervals(rbind(c(-Inf, 0), c(-C/B, Inf))))
            }
          } else {
            warning("Ill-conditioned quadratic")
            return(Intervals(c(-Inf,0)))
          }
        }
      }
    })
    }
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

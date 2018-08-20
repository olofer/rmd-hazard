#
# ASSORTED SUBPROGRAMS FOR SDE & PDE ILLUSTRATIONS
#
#

#
# The below 2 versions of value sampling presumably provide the same mean values
# but they have different variances;
#

sim.sde.euler.kill <- function(sdep, t0, x0, T, L, dt, save.path = FALSE) {
  t <- t0
  x <- x0
  uw <- sqrt(dt)
  maxsteps <- round((T - t0) / dt)
  if (save.path) {
    tvec <- array(NA, maxsteps)
    xvec <- array(NA, maxsteps)
    fvec <- array(NA, maxsteps)
  }
  is.absorbed <- (abs(x) >= L)
  is.killed <- FALSE
  stopifnot(!is.absorbed)
  utot <- 0
  jj <- 0
  while (jj < maxsteps && !is.absorbed && !is.killed) {
    dw <- uw * rnorm(1)
    dx <- sdep$mu(x, t) * dt + sdep$sigma(x, t) * dw
    f <- sdep$f(x, t)
    jj <- jj + 1
    if (save.path) {
      tvec[jj] <- t
      xvec[jj] <- x
      fvec[jj] <- f
    }
    utot <- utot + dt * f
    z <- exp(-dt * sdep$h(x, t))  # probability of survival for this interval
    x <- x + dx
    t <- t + dt
    is.absorbed <- (abs(x) >= L)
    is.killed <- (runif(1) >= z)
  }
  if (!is.absorbed && !is.killed) {
    utot <- utot + sdep$psi(x)
  }
  if (save.path) {
    list(
      t = tvec[1:jj],
      x = xvec[1:jj],
      f = fvec[1:jj],
      utot = utot,
      nsteps = jj,
      is.absorbed = is.absorbed,
      is.killed = is.killed)
  } else {
    list(
      t = t,
      x = x,
      utot = utot,
      nsteps = jj,
      is.absorbed = is.absorbed,
      is.killed = is.killed)
  }
}

sim.sde.euler <- function(sdep, t0, x0, T, L, dt, save.path = FALSE) {
  t <- t0
  x <- x0
  maxsteps <- round((T - t0) / dt)
  uw <- sqrt(dt)
  if (save.path) {
    tvec <- array(NA, maxsteps)
    xvec <- array(NA, maxsteps)
    fvec <- array(NA, maxsteps)
    hvec <- array(NA, maxsteps)
    pvec <- array(NA, maxsteps)
    uvec <- array(NA, maxsteps)
  }
  is.absorbed <- (abs(x) >= L)
  stopifnot(!is.absorbed)  # disallow starting outside domain
  hcum <- 0
  ucum <- 0
  jj <- 0
  while(jj < maxsteps && !is.absorbed) {
    jj <- jj + 1
    f <- sdep$f(x, t)
    h <- sdep$h(x, t)
    p <- exp(-hcum)
    if (save.path) {
      tvec[jj] <- t
      xvec[jj] <- x
      fvec[jj] <- f
      hvec[jj] <- h
      pvec[jj] <- p
      uvec[jj] <- ucum
    }
    is.absorbed <- (abs(x) >= L)
    if (is.absorbed) break
    dw <- uw * rnorm(1)
    dx <- sdep$mu(x, t) * dt + sdep$sigma(x, t) * dw
    x <- x + dx
    hcum <- hcum + h * dt
    ucum <- ucum + p * f * dt
    t <- t + dt
  }
  if (!is.absorbed) ucum <- ucum + p * sdep$psi(x)
  if (save.path) {
    list(
      t = tvec[1:jj],
      x = xvec[1:jj],
      f = fvec[1:jj],
      h = hvec[1:jj],
      p = pvec[1:jj],
      u = uvec[1:jj],
      utot = ucum,
      nsteps = jj,
      is.absorbed = is.absorbed)
  } else {
    list(
      t = t,
      x = x,
      utot = ucum,
      nsteps = jj,
      is.absorbed = is.absorbed)
  }
}

#
# Below follows a few routines for sample-path post-processing "accounting"
#

# Initialize a counting "object" with nx-by-nt grid cells
create.xt.cells <- function(nx, x12, nt, t12) {
  xbrk <- seq(from = x12[1], to = x12[2], len = nx + 1)
  xmid <- (xbrk[2:(nx + 1)] + xbrk[1:nx]) / 2
  tbrk <- seq(from = t12[1], to = t12[2], len = nt + 1)
  tmid <- (tbrk[2:(nt + 1)] + tbrk[1:nt]) / 2
  list(
    xbrk = xbrk, xmid = xmid,
    tbrk = tbrk, tmid = tmid,
    N = array(0, c(nx, nt)),
    U = array(0, c(nx, nt)),
    W = array(-1, c(nx, nt))
    )
}

# Update the counter object with one path sample
update.xt.cell <- function(xto, r) {
  nt <- length(xto$tmid)
  nx <- length(xto$xmid)
  nk <- length(r$t)
  idxh <- match('h', names(r))
  use.explicit.hazard <- !is.na(idxh)
  dt <- r$t[2] - r$t[1]
  kkp <- 1
  zk <- r$utot
  for (jj in 1:nt) {
    if (xto$tbrk[jj + 1] <= r$t[1] || xto$tbrk[jj] > r$t[nk]) next
    # find the time-index kk in r$t that is closest to the query time tmid[jj]
    kk <- match(TRUE, r$t >= xto$tmid[jj])
    #stopifnot(!is.na(kk))
    if (is.na(kk)) next
    stopifnot(kk >= kkp)
    stopifnot(r$t[kk] >= xto$tbrk[jj] && r$t[kk] < xto$tbrk[jj + 1])
    # then locate the space bin ii that r$x[kk] falls in
    ii <- match(TRUE, xto$xbrk[1:nx] <= r$x[kk] & xto$xbrk[2:(nx + 1)] > r$x[kk])
    if (is.na(ii)) next
    # then evaluate the future integral sample vij by "unrolling the total integral"
    # z <- r$utot; for (tt in 1:kk) z <- (z - dt * ft) * exp(dt * ht); vij <- z
    if (use.explicit.hazard) {
      while (kkp < kk) {
        zk <- (zk - dt * r$f[kkp]) * exp(dt * r$h[kkp])
        kkp <- kkp + 1
      }
    } else {
      while (kkp < kk) {
        zk <- (zk - dt * r$f[kkp])
        kkp <- kkp + 1
      }
    }
    vij <- zk
    nij <- xto$N[ii, jj]
    uij <- xto$U[ii, jj]
    wij <- xto$W[ii, jj]
    bet <- nij / (nij + 1)
    uij <- uij * bet + vij * (1 - bet)
    wij <- wij * bet + vij * vij * (1 - bet)
    nij <- nij + 1
    # store {uij, wij, nij} back into xto
    xto$N[ii, jj] <- nij
    xto$U[ii, jj] <- uij
    xto$W[ii, jj] <- wij
    kkp <- kk
  }
  xto
}

# Merge two counter objects; error out if incompatible
merge.xt.cell <- function(xto1, xto2) {
  stopifnot(length(xto1$xbrk) == length(xto2$xbrk))
  stopifnot(all(xto1$xbrk == xto2$xbrk))
  stopifnot(length(xto1$tbrk) == length(xto2$tbrk))
  stopifnot(all(xto1$tbrk == xto2$tbrk))
  xto <- xto1
  N12 <- xto1$N + xto2$N
  A12 <- xto1$U * xto1$N + xto2$U * xto2$N
  B12 <- xto1$W * xto1$N + xto2$W * xto2$N
  xto$N <- N12
  nodata12 <- which(N12 == 0)
  A12 <- A12 / N12
  A12[nodata12] <- 0
  xto$U <- A12
  B12 <- B12 / N12
  B12[nodata12] <- -1
  xto$W <- B12
  xto
}

# Solve the backward equation for u(t, x) over t = T...t0 with terminal
# condition u(T, x) = psi(x) and with  absorbing boundary u(t, L) = u(t,-L) = 0
# n is the Chebyshev grid size; dt is the time-step (implicit Euler method)
solve.backward.pde <- function(sdep, t0, T, L, dt, n) {
  stopifnot(is.list(sdep))
  stopifnot(t0 < T)
  c <- cheby.cheb(n)
  ell <- L - (-L)
  y <- -L + (ell / 2) * (c$x + 1) # y in [-L, L] from x in [-1, 1]
  Dy <- (2 / ell) * c$D
  Dyy = Dy %*% Dy
  nk = round((T - t0) / dt)
  U <- array(NA, c(n - 1, nk))  # buffer for collocation representation of solution
  tU <- array(NA, nk)
  kk <- nk
  tk <- T
  uk <- sdep$psi(y)
  uk <- uk[2:n] # end points are enforced hard implicit zeros 
  while (kk >= 1) {
    U[, kk] <- uk
    tU[kk] <- tk
    Ak <- diag(sdep$mu(y, tk)) %*% Dy + diag(0.5 * sdep$sigma(y, tk) ^ 2) %*% Dyy - diag(sdep$h(y, tk))
    Ak <- Ak[2:n, 2:n] # due to the absorbing b.c.
    fk <- sdep$f(y, tk)
    qr.of.Ak <- qr(diag(n - 1) - dt * Ak)
    uk <- qr.solve(qr.of.Ak, uk + dt * fk[2:n])
    kk <- kk - 1
    tk <- tk - dt
  }
  list(t = tU, U = U, x = y, L = L)
}

# Solve the forward (Fokker-Planck) equation for t = t0..T
# with the probability distribution for X at t=t0 given by a Gaussian
# centered at x0 and with standard deviation sigma0:
# p_t = -(mu*p)_x + 0.5*(sigma^2*p)_xx - h*p; p(x,t=t0) ~ N(x-x0, sigma0^2)
solve.forward.pde <- function(sdep, t0, T, L, dt, n, x0, sigma0, int.u0 = FALSE) {
  stopifnot(is.list(sdep))
  stopifnot(t0 < T)
  c <- cheby.cheb(n)
  ell <- L - (-L)
  y <- -L + (ell / 2) * (c$x + 1)
  Dy <- (2 / ell) * c$D
  Dyy = Dy %*% Dy
  if (int.u0) {
    cc <- cheby.clencurt(n)
    w <- (ell / 2) * cc$w
  }
  nk = round((T - t0) / dt)
  P <- array(NA, c(n - 1, nk))
  tP <- array(NA, nk)
  kk <- 0
  tk <- t0
  pk <- exp(-(y - x0)^2 / (2*sigma0^2)) / (sqrt(2*pi)*sigma0)
  pk <- pk[2:n]
  u0 <- 0
  while (kk < nk) {
    kk <- kk + 1
    P[, kk] <- pk
    tP[kk] <- tk
    if (int.u0) {
      u0 <- u0 + dt * as.vector(w %*% (c(0, pk, 0) * sdep$f(y, tk)))  # integrate expected value at this time
    }
    Ak <- -Dy %*% diag(sdep$mu(y, tk)) + 0.5 * Dyy %*% diag(sdep$sigma(y, tk) ^ 2) - diag(sdep$h(y, tk))
    Ak <- Ak[2:n, 2:n]
    qr.of.Ak <- qr(diag(n - 1) - dt * Ak)
    pk <- qr.solve(qr.of.Ak, pk)
    tk <- tk + dt
  }
  if (int.u0) {
    u0 <- u0 + as.vector(w %*% (c(0, pk, 0) * sdep$psi(y)))  # add terminal value
  }
  list(t = tP, U = P, x = y, L = L, u0 = u0)
}

# Post-process the returned struct from the previous function(s)
render.pde <- function(pdep, xx) {
  stopifnot(is.list(pdep))
  nt <- length(pdep$t)
  nx <- length(xx)
  Uxx <- array(0, c(nx, nt))
  ell <- 2 * pdep$L  # assume range [-L, L]
  yy <- 2 * (xx + pdep$L) / ell - 1
  stopifnot(all(yy >= -1 & yy <= 1))
  for (kk in 1:nt) {
    Uxx[, kk] <- cheby.bary(c(0, pdep$U[, kk], 0), yy)
  }
  list(x = xx, t = pdep$t, Uxt = Uxx)
}

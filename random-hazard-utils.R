#
# Two types of utility functions:
#   1. Chebyshev polynomial tools (first in file)
#   2. Some first-visit counting tools (further down in file)
#
#
# Author: K Erik J Olofsson
#


#
# Utility functions for Chebyshev polynomial approximation.
# Differentiation, interpolation, and quadrature.
# Basic usages is for the bounded interval [-1, +1].
# Unbounded domains can be transformed as in Boyd's book (Dover, 2001).
# Note that Trefethen's functions list nodes x in decreasing order.
#

#
# Chebyshev nodes x and the associated differentiation matrix D
#
cheby.cheb <- function(n) {
  if (n == 0) {
    list(x = 1, D = 0)
	} else {
    x <- cos(pi*(0:n)/n)
    c <- c(2, array(1, n - 1), 2) * (-1)^(0:n)
    X <- array(x, c(n + 1, n + 1))
    dX <- X - t(X)
    D <- outer(c, 1 / c) / (dX + diag(1, n + 1))
    D <- D - diag(rowSums(D), n + 1)
    list(x = x, D = D)
  }
}

#
# Quadrature weights w for nodes x
#
cheby.clencurt <- function(n) {
  theta <- pi * (0:n) / n
  x <- cos(theta)
  w <- array(0, n + 1)
  ii <- 2:n
  v <- array(1, n - 1)
	if (n %% 2 == 0) {
    w[1] <- 1 / (n^2 - 1)
    w[n + 1] <- w[1]
    for (kk in 1:(n / 2 - 1)) {
      v <- v - 2 * cos(2 * kk * theta[ii]) / (4 * kk^2 - 1)
    }
    v <- v - cos(n * theta[ii]) / (n^2 - 1)
  } else {
    w[1] <- 1 / n^2
    w[n + 1] <- w[1]
    for (kk in 1:((n - 1) / 2)) {
      v <- v - 2 * cos(2 * kk * theta[ii])/(4 * kk^2 - 1)
    }
  }
  w[ii] <- 2 * v / n
  list(x = x, w = w)
}

#
# Barycentric interpolation; see Berrut, SIAM Review 2004
# xx should be in range [-1, +1]
#
cheby.bary <- function(f, xx) {
  n <- length(f) - 1
  x <- cos(pi * (0:n) / n)
  cc <- c(1/2, array(1, n - 1), 1/2) * (-1)^(0:n)
  nn <- length(xx)
  numer <- array(0, nn)
  denom <- array(0, nn)
  exact <- array(0, nn)
  for (ii in 1:(n + 1)) {
    xdiff <- xx - x[ii]
    tmp <- cc[ii] / xdiff
    numer <- numer + tmp * f[ii]
    denom <- denom + tmp
    exact[xdiff == 0] <- ii 
  }
  ff <- numer / denom
  jj <- which(exact > 0)
  if (length(jj) > 0) {
    ff[jj] <- f[exact[jj]]
  }
  ff  ## interpolated value of (x, f) in grid xx
}

#
# Solve the PDE:
#
#   u_t = a(x) * u_x + 0.5 * b(x)^2 * u_xx - h(x) * u
#
# with initial data u(t = 0, x) = 1 on the interval [xa, xb]
# with absorbing boundary conditions u(t, xa) = u(t, xb) = 0, 
# for 0 <= t <= tmax.
#
# Use Chebyshev node space discretization to obtain u_t = A * u
# This is solved in time with backward Euler stepping:
#
# u_1 - u_0 = dt * A * u_1, so (I - dt * A) * u_1 = u_0
#
cheb.solve.PDE.absorbing <- function(ax, bx, hx, n, dt, xa, xb, tmax, xx = NULL) {
	stopifnot(is.function(ax) && is.function(bx) && is.function(hx))
	stopifnot(xa < xb && dt > 0 && tmax > dt)
	c <- cheby.cheb(n)
	ell <- xb - xa
	y <- xa + (ell / 2) * (c$x + 1) # y in [xa, xb] from x in [-1, 1]
	Dy <- (2 / ell) * c$D
	Dyy = Dy %*% Dy
	A <- diag(ax(y)) %*% Dy + diag(0.5 * bx(y) ^ 2) %*% Dyy - diag(hx(y))
	A <- A[2:n, 2:n] # due to the absorbing b.c.
	# Setup backward Euler propagator; factorize dense matrix I - dt * A
	qr.of.A <- qr(diag(n - 1) - dt * A)
	# Initialize 
	k <- 0; t <- 0; u <- array(1, n - 1)
	nmax <- round(tmax / dt)
	T <- array(0, nmax)
	U <- array(0, c(nmax, n - 1))
	UX <- NA
	if (!is.null(xx)) {
    UX <- array(0, c(nmax, length(xx)))
    yy <- 2 * (xx - xa) / ell - 1
    stopifnot(all(yy >= -1 & yy <= 1))
  }
	while (k < nmax) {
		k <- k + 1
		T[k] <- t
		U[k, ] <- u
		if (!is.null(xx)) {
			UX[k, ] <- cheby.bary(c(0, u, 0), yy)
		}
		u <- qr.solve(qr.of.A, u)
		t <- t + dt
	}
	list(t = T, uc = U, ux = UX)
}

#
# Solve the PDE:
#
#   u_t = a(x) * u_x + 0.5 * b(x)^2 * u_xx - h(x) * u
#
# with initial data u(t = 0, x) = 1 on the entire real line x, 
# for 0 <= t <= tmax.
#
# This is done with the rational Chebyshev domain transformation.
# See Boyd's book on Pseudospectral methods for reference.
#
# L is a length scale parameter for the domain transformation:
#
#   y = x / sqrt(L^2 + x^2)
#
# The new coordinate y covers the domain [-1, 1].
#
# REQUIRE: hx(-Inf) = finite AND hx(Inf) = finite
# (it is assumed that u(t, x) should decay to zero as |x| goes large)
#
cheb.solve.PDE.unbounded <- function(ax, bx, hx, n, dt, L, tmax, xx = NULL) {
  stopifnot(is.function(ax) && is.function(bx) && is.function(hx))
  stopifnot(L > 0 && dt > 0 && tmax > dt)
  c <- cheby.cheb(n)
  Q <- 1 - c$x^2
  y <- L * c$x / sqrt(Q)  # NOTE: y[1] = -Inf and y[n+1] = Inf here
  tmp1 <- hx(y)
  stopifnot(all(is.finite(tmp1)))
  A1 <- diag(-tmp1)  # -h*u
  tmp2 <- ax(y); tmp2[1] <- 0; tmp2[n + 1] <- 0  # need to make sure Inf's are not retained
  stopifnot(all(is.finite(tmp2)))
  tmp3 <- bx(y); # tmp3[1] <- 0; tmp3[n + 1] <- 0  # only adjust if required
  stopifnot(all(is.finite(tmp3)))
  A2 <- (diag( (1/L) * sqrt(Q) * Q * tmp2) - diag((3/2) * (1/L^2) * tmp3^2 * Q^2 * c$x)) %*% c$D
  A3 <- diag((1/2) * (1/L^2) * tmp3^2 * Q^3) %*% c$D %*% c$D
  A <- A1 + A2 + A3
  qr.of.A <- qr(diag(n + 1) - dt * A)
#  M <- solve(diag(n + 1) - dt * A)
  nt <- round(tmax / dt)
  U <- array(0, c(nt, n + 1))
  T <- array(0, nt)
  UX <- NA
  if (!is.null(xx)) {
    UX <- array(0, c(nt, length(xx)))
    yy <- xx / sqrt(L^2 + xx^2)
    stopifnot(all(is.finite(yy)))
  }
  u <- array(1, n + 1)
  k <- 0; t <- 0
  while (k < nt) {
    k <- k + 1
    T[k] <- t
    U[k, ] <- u
    if (!is.null(xx)) {
      UX[k, ] <- cheby.bary(u, yy)
    }
    u <- qr.solve(qr.of.A, u)
  #  u <- M %*% u
    t <- t + dt
  }
  list(t = T, uc = U, ux = UX)
}

#
# SURPRISE BONUS CODE
#
# Test function to illustrate how to solve the 1D quantum
# harmonic oscillator on the *infinite* real line using variable transformation
# plot the first few eigenfunctions, normalized using Clenshaw-Curtis quadrature.
# Return the eigenvalues: should be 1/2, 3/2, 5/2 etc..
# Plot an an equidistant grid in y, using barycentric interpolation in x
# example: E <- cheb.qho.1d(60, 3, 8, 4)
#
# The trick below (and above) is to use the following transformations.
# y is [-inf,inf], x is [-1,+1]
# Q = 1-x^2, y=L*x/sqrt(Q), x=y/sqrt(L^2+y^2), dy=L*dx/Q^(3/2)
# u_y = sqrt(Q)*Q*u_x/L
# u_yy = Q^2*{Q*u_xx-3*x*u_x}/L^2
# also: y=L*cot(t), t=arccot(y/L), t=cos(x), t=arccos(x)
# smallest delta-x is L/n
#
cheb.qho.1d <- function(n, L, neig, Ymax) {
  C <- cheby.cheb(n)
  W <- cheby.clencurt(n)
  stopifnot(L > 0)
  stopifnot(n > 0)
  stopifnot(Ymax > 0)
  stopifnot(neig >= 0)
  D2 <- C$D %*% C$D
  Q <- 1 - C$x^2
  w <- L*(W$w) / (sqrt(Q)^3) # yep, do not forget to transform integral weight too
  A <- (-(1/(2*L^2)) * diag(Q^3, n + 1) %*% D2
    + (3/(2*L^2)) * diag(Q^2 * C$x, n + 1) %*% C$D
    + (L^2/2) * diag(C$x^2/Q, n + 1) )
  A <- A[2:n, 2:n]	# implied b.c. is u->0 at both infinities
  E <- eigen(A)
  colvec <- rainbow(neig + 1)
  yy <- seq(from = -Ymax, to = Ymax, by = 2*Ymax/1000)
  xx <- yy / sqrt(yy^2+L^2)
  eigk <- c(0, E$vectors[, n-1], 0)
  eigk <- eigk / sqrt(sum(w[2:n] * abs(eigk[2:n])^2))
  ff <- cheby.bary(eigk, xx)
  fmax <- 1 / pi^(1/4) # should be the max value for the 0-th eigenfcn; check this
  plot(x = yy, y = ff, type = "l",
    col = colvec[neig + 1], xlab = "y", ylab = "QHO eigenfunctions",
    main = sprintf("First %i 1D QHO eigenfunctions (Hermite fcns)", neig + 1),
    xlim = c(-Ymax, Ymax), ylim=c(-fmax, fmax))
  for (jj in 1:neig) {
    eigk <- c(0, E$vectors[, n-1-jj], 0)
    eigk <- eigk / sqrt(sum(w[2:n] * abs(eigk[2:n])^2))
    ff <- cheby.bary(eigk, xx)
    lines(x = yy, y = ff, type = "l", col = colvec[jj])
  }
  lines(x = c(-Ymax, Ymax), y = c(fmax, fmax), col = "black")
  lines(x = c(-Ymax, Ymax), y = c(-fmax, -fmax), col = "black")
  E$values
}

#
# Aux. functions for counting statistics.
# These are designed to work with foreach loop constructs.
#

#
# (t, x) is a path, e is logical TRUE if (t, x) ends with an event, otherwise FALSE
# xg is a monotonic sequence defining a grid for covariate x.
# tg is a monotonic sequence defining a grid for time t.
#
get.count.matrices <- function(t, x, e, xg, tg) {
  nt <- length(t)
  stopifnot(!is.unsorted(t))
  stopifnot(length(x) == nt)
  stopifnot(length(e) == 1 && is.logical(e))
  nx <- length(xg) - 1
  stopifnot(nx > 0)
  stopifnot(!is.unsorted(xg))
  ntg <- length(tg)
  stopifnot(tg[1] == 0)
  stopifnot(!is.unsorted(tg))
  cnt <- array(0, c(nx, ntg))
  cnd <- array(0, c(nx, ntg))
  for (i in 1:nx) {
    idx <- which(x >= xg[i] & x <= xg[i + 1])
    if (length(idx) == 0) next
    idx.first <- min(idx)
    ti <- t[idx.first:nt] - t[idx.first]  # future time going forward
    ti.end <- ti[length(ti)]
    for (j in 1:ntg) {
      if (tg[j] <= ti.end) {
        cnt[i, j] <- cnt[i, j] + 1
      } else {
        if (e) {
          cnd[i, j] <- cnd[i, j] + 1
        }
      }
    }
  }
  list(cnt = cnt, cnd = cnd)
}


# .multicombine = FALSE
combine.count.matrices <- function(a, b) {
  list(cnt = a$cnt + b$cnt, cnd = a$cnd + b$cnd)
}

# .multicombine = TRUE
multicombine.count.matrices <- function(...) {
  alist <- list(...)
  stopifnot(length(alist) > 0)
  mt <- alist[[1]]$cnt
  md <- alist[[1]]$cnd
  for (j in 2:length(alist)) {
    mt <- mt + alist[[j]]$cnt
    md <- md + alist[[j]]$cnd
  }
  list(cnt = mt, cnd = md)
}

extract.survival.map <- function(a, novalue = 0) {
  nx <- nrow(a$cnt)
  nt <- ncol(a$cnt)
  srv <- array(NA, c(nx, nt))
  srv[, 1] <- 1
  idx.empty <- which(a$cnt[, 1] == 0)
  if (length(idx.empty) != 0) {
    srv[idx.empty, 1] <- novalue
  }
  for (j in 2:nt) {
    dj <- a$cnd[, j] - a$cnd[, j - 1]  # num. events in last interval
    nj <- a$cnt[, j - 1]  # num. "at risk" in last interval
    stopifnot(all(dj <= nj))  # sanity check
    pvj <- array(1, c(nx, 1)) - dj / nj
    if (any(nj == 0)) {
      pvj[which(nj == 0), 1] <- 1
    }
    srv[, j] <- pvj * srv[, j - 1]
  }
  srv
}

get.midpoints <- function(xg) {
  stopifnot(!is.unsorted(xg))
  nx <- length(xg) - 1
  xvec <- array(NA, c(nx, 1))
  for (j in 1:nx) {
    xvec[j] <- (xg[j] + xg[j + 1]) / 2
  }
  xvec
}

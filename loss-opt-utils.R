#
# Subprograms for "generalized linear model" type loss minimization.
# Uses second order optimization (Newton-type) algorithm.
#
# The loss function is customizable but should be convex.
# L2 regularization is supported.
#
# A useful feature is warm-started re-optimization for scans
# over a progression of L2 parameters.
#
#
# K Erik J Olofsson, February/March 2018
#

lghfunc.ell2 <- function(y, f, w = NULL) {
  e <- y - f
  if (is.null(w)) {
    l <- 0.5 * e ^ 2
    g <- e * (-1)
    h <- array(1, length(y))
  } else {
    l <- 0.5 * w * e ^ 2
    g <- e * (-w)
    h <- w
  }
  list(loss = l, grad = g, hess = h)
}

lghfunc.xentropy <- function(y, f, w = NULL) {
  stopifnot(is.null(w))
  z <- 1 / (1 + exp(-f))
  l <- -(y * log(z) + (1 - y) * log(1 - z))
  g <- z - y
  h <- z * (1 - z)
  list(loss = l, grad = g, hess = h)
}

# NOTE: xentlambda degenerates to xentropy for unit weights w (or NULL)
lghfunc.xentlambda <- function(y, f, w = NULL) {
  if (is.null(w)) {
    return(lghfunc.xentropy(y, f, w))
  } else {
    z <- 1 - (1 + exp(f))^(-w)
    l <- -(1 - y) * log(1 - z) - y * log(z)
    xi <- 1 / (1 + exp(-f))
    g <- xi * w * (1 - y / z)
    x <- 1 + exp(f)
    ff <- w * exp(-f) *xi ^ 2;
    gg <- (x^w * (1 + w * exp(f) - x^w)) / ((x^w - 1)^2)
    h <- ff * (1 + y * gg)
  }
  list(loss = l, grad = g, hess = h)
}

#
# Fit function on the form
#
#   fhat = X %*% beta
#
# by optimizing parameter vector beta; where the total objective is:
#
#   J(beta) = (1/n) * sum_{i=1..n} ow[i] * lghfunc(y, b + fhat, w)$loss[i] 
#             + (1/2) * sum_{j=1..m} ell2[j] * beta[j]^2
#
# This is solved by a (second order) Newton iteration method.
#
# lghfunc is the pointwise loss (also returns gradient grad and hessian hess)
# y are labels / targets, one per example
# X is the n-by-m matrix of covariates; one example per row
# w are optional example weights (generic, passed onto custum loss function)
#
# b is an optional (fixed) bias vector (useful when boosting with adapting basis)
# ow is an optional linear observation weight vector
# (note: linear weights can also be specified with w, depending on lghfunc)
#
# ell2 is a vector of regularization coefficients (one per column of X)
# (scalar ell2 is auto-expanded to m-vector)
# (ell2 can also be a matrix of size m-by-m; reg = 0.5 * t(beta) %*% ell2 %*% beta)
#
# beta.init is the initial coefficient guess for warm-starting (zero by default)
# eta is a Newton method damping scalar
# eptol is a convergence requirement / stop condition
# maxiter is the maximum number of allowed Newton iterations
# beta.iter is a boolean: store the iteration progress for beta
#
# stop.grad and stop.hess are booleans indicating whether gradient
# or hessian at the final iteration should be returned (FALSE default)
#
opt.loss.ell2 <- function(
  lghfunc,
  y, 
  X, 
  w = NULL, 
  b = NULL,
  ow = NULL, 
  ell2 = 0,
  beta.init = array(0, ncol(X)), 
  eta = 1.0, 
  eptol = 1e-10, 
  maxiter = 100,
  beta.iter = FALSE,
  stop.grad = FALSE,
  stop.hess = FALSE)  # return the "data hessian" (not the regularized loss hessian)
{
  # First do some type-checking / documentation of intent
  stopifnot(is.function(lghfunc))
  stopifnot(is.matrix(X))
  n <- nrow(X)
  m <- ncol(X)
  stopifnot(length(y) == n)
  if (!is.null(w)) {
    stopifnot(length(w) == n)
  }
  if (!is.null(b)) {
    if (length(b) == 1) b <- array(b, n)  # auto-expand scalar
    stopifnot(length(b) == n)
  }
  if (!is.null(ow)) {
    stopifnot(length(ow) == n)
    stopifnot(all(ow >= 0))
    sow <- sum(ow)
    stopifnot(sow > 0)
    if (sow != 1) {  # renormalize
      ow <- ow / sow
    }
  }

  # Further sanity checks
  stopifnot(eta > 0 && eta <= 1)
  stopifnot(eptol >= 0)
  stopifnot(maxiter > 0)

  if (length(ell2) == 1 && m > 1) ell2 <- array(ell2, m)  # auto-expand scalar, if needed

  if (is.matrix(ell2)) {
    stopifnot(ncol(ell2) == m && nrow(ell2) == m)
  } else {
    stopifnot(length(ell2) == m)
    stopifnot(all(ell2 >= 0))
  }

  stopifnot(length(beta.init) == m)

  if (beta.iter) {
    beta.log <- array(NA, c(maxiter, m))
  }

  beta <- beta.init
  kk <- 0

  # allocate storage for progress information
  obj.tot <- array(NA, maxiter)
  obj.fit <- array(NA, maxiter)
  gnorm <- array(NA, maxiter)

  lnorm0 = max(1, max(abs(ell2)))  # max(abs(.)) is the inf-norm

  while (TRUE) {
    fhat <- as.vector(X %*% beta)
    if (!is.null(b)) fhat <- fhat + b
    lgh <- lghfunc(y, fhat, w)
    if (!is.null(ow)) {
      lgh$loss <- lgh$loss * ow
      lgh$grad <- lgh$grad * ow
      lgh$hess <- lgh$hess * ow
    }
    Jfit <- (1/n) * sum(lgh$loss)
    if (is.matrix(ell2)) {
      Jreg <- (1/2) * (beta %*% (ell2 %*% beta))
      gvec <- (1/n) * apply(X, 2, `%*%`, lgh$grad) + ell2 %*% beta
    } else {
      Jreg <- (1/2) * (ell2 %*% (beta ^ 2))
      gvec <- (1/n) * apply(X, 2, `%*%`, lgh$grad) + ell2 * beta
    }
    stopifnot(length(gvec) == m)
    obj.fit[kk + 1] <- Jfit 
    obj.tot[kk + 1] <- Jfit + Jreg
    gnorm[kk + 1] <- max(abs(gvec))
    if (beta.iter) {
      beta.log[kk + 1, ] <- beta
    }
    is.converged <- (gnorm[kk + 1] / lnorm0) < eptol
    is.done <- (kk == maxiter) || is.converged
    kk <- kk + 1
    if (is.done && !stop.hess) break  # if hessian is to be returned proceed ...
    Z <- apply(X, 2, `*`, lgh$hess)
    if (is.done) break # .. then stop
    if (is.matrix(ell2) || length(ell2) == 1) {
      H <- (1/n) * (t(Z) %*% X) + ell2
    } else {
      H <- (1/n) * (t(Z) %*% X) + diag(ell2)
    }
    stopifnot(nrow(H) == m && ncol(H) == m)
    # TODO: switchable methods: Cholesky, solve, QR; or explicitly the IRWLS "lookalike"
    beta <- beta - eta * solve(H, as.vector(gvec))  # take (damped) Newton step
  }

  list(
    beta = beta,
    is.converged = is.converged,
    iters = kk,
    obj = obj.tot[1:kk],
    obj.fit = obj.fit[1:kk],
    gnorm = gnorm[1:kk],
    beta.log = if(beta.iter) beta.log[1:kk, ] else NA,
    grad = if (stop.grad) gvec else NA,
    hess = if(stop.hess) { t(Z) %*% X } else { NA }
    )
}

#
# Univariate "smooth curve" fitting application of the above "glm" type sub-program.
# Covariate x need not be sorted.
#
# Fits a piecewise constant "curve" where each delta is regularized according to
# either 0-th, 1st-, or 2nd-order ell-2 penalty.
#
# regtyp = 0, 1, or 2 (ignored if ell2.vec is NULL)
# 
# It is also possible to specify a "prior" observation with yreg.
#
# Uncertainty region is automatically generated from the diagonal of the Hessian
# matrix at the fitted parameter values (Fisher Information)
#
opt.loss.ell2.pwco.1d <- function(
  y,
  x,
  w = NULL,
  nn = 30,
  xrange = c(min(x), max(x)),
  ell2.vec = c(0),
  regtyp = 2,
  yreg = NULL,
  lghfunc = lghfunc.xentlambda,
  warm.restart = FALSE)
{
  n <- length(y)
  stopifnot(length(x) == n)
  if (!is.null(w)) {
    stopifnot(length(w) == n)
  }
  stopifnot(length(xrange) == 2)
  stopifnot(xrange[2] > xrange[1])
  xg <- seq(from = xrange[1], to = xrange[2], len = nn + 1)  # mini-interval edges (defines the basis set)
  X <- array(NA, c(n, nn))  # the covariate vector (one basis per column)
  for (ii in 1:nn) {
    X[, ii] <- ifelse(x > xg[ii] & x <= xg[ii + 1], 1, 0)
  }
  if (!is.null(yreg)) {
    # TODO: might need extended & improved options for setting the interval based prior observations
    # "prior observations" / useful to maintain reasonable estimates even when there are very few y = 1 labels
    stopifnot(length(yreg) == 1)
    if (yreg >= 0 && yreg <= 1) {
      # here insert given y = yreg, mean(w) observations
      yadd <- array(yreg, nn)
    } else {
      # special auto-tune based on "average intensity"
      yadd <- array(mean(y), nn)
    }
    Xadd <- diag(array(1, nn))
    X <- rbind(X, Xadd)
    y <- c(y, yadd)
    if (!is.null(w)) {
      wadd <- array(mean(w), nn)
      w <- c(w, wadd)
    }
#    warning('yreg argument not handled')
  }
  nell <- length(ell2.vec)
  stopifnot(nell > 0)
  stopifnot(all(ell2.vec >= 0))
  stopifnot(regtyp == 0 || regtyp == 1 || regtyp == 2) 
  if (regtyp == 1) {
    D.op <- array(0, c(nn - 1, nn))
    for (ii in 1:nrow(D.op)) {
      D.op[ii, ii:(ii+1)] <- c(-1, 1);
    }
  } else if (regtyp == 2) {
    D.op <- array(0, c(nn - 2, nn))
    for (ii in 1:nrow(D.op)) {
      D.op[ii, ii:(ii+2)] <- c(-1, 2, -1);
    }
  } else {
    D.op <- diag(array(1, nn))
  }
  L <- t(D.op) %*% D.op  # reg. matrix to be scaled by lambda >= 0
  opt.list <- list() 
  for (ii in 1:nell) {
    # here estimate one model per lambda, with optional warm-starting
    # NOTE: warm start may only make sense if ell2.vec is sorted (but this is NOT checked for)
    if (ii > 1 && warm.restart && opt.list[[ii - 1]]$is.converged) {
      beta.init <- opt.list[[ii - 1]]$beta
    } else {
      # use zero vector as initial beta vector
      beta.init <- array(0, nn)
    }
    # single instance solution
    opt.list[[ii]] <- opt.loss.ell2(
      lghfunc, y, X, w, 
      b = NULL, ow = NULL, 
      ell2 = as.matrix(ell2.vec[ii] * L),
      beta.init = beta.init,
      eta = 1.0, eptol = 1e-10, maxiter = 100,
      beta.iter = FALSE, stop.grad = TRUE, stop.hess = TRUE)
    if (!opt.list[[ii]]$is.converged) {
      warning(sprintf('not converged for ell2 = %e', ell2.vec[ii]))
    }
  }
  # Return list of return lists (1 per ell2.vec element)
  opt.list
}

#
# K-fold CV tuner for selection of reasonable ell2 parameter value given a dataset.
# It will be parallelized if the foreach environment is specified (outside of this function)
#
auto.pwco.1d <- function(y, x, w = NULL, xrange = NULL) {
  return(NA)
}

#
# Plot the family of solutions provided by the function
# "opt.loss.ell2.pwco.1d" above.
#
plot.pwco.1d <- function(
  ol,
  xrange = NULL,
  ytyp = 0,
  dt = NULL,
  FisherInfo = FALSE) 
{
  # ol is the returned list of lists from "opt.loss.ell2.pwco.1d" above
  stopifnot(is.list(ol))
  nsols <- length(ol)
  stopifnot(nsols >= 1)
  m <- length(ol[[1]]$beta)
  if (is.null(xrange)) {
    xvec <- 1:m
    xlab <- 'coef index'
  } else {
    stopifnot(length(xrange) == 2)
    xseq <- seq(from = xrange[1], to = xrange[2], len = m + 1)  # use 'midpoints'
    xvec <- array(NA, m)
    for (ii in 1:m) {
      xvec[ii] <- sum(xseq[ii:(ii + 1)]) / 2
    }
    xlab <- 'x'
  }
  if (is.null(dt)) dt <- 1
  stopifnot(dt >= 0)
  stopifnot(ytyp == 0 || ytyp == 1)
  if (ytyp == 0) {
    ylab <- 'log-odds f (beta)'
    yfunc <- function(a) { a }
  } else if (ytyp == 1) {
    ylab <- 'hazard h'
    yfunc <- function(a) { log(1 + exp(a)) / dt }
  }
  lwd <- 2
  col1 <- 'blue'
  if (nsols > 1) {
    require(viridis)
    coll <- viridis(nsols)
    col1 <- coll[1]
  }
  plot(
    x = xvec, y = yfunc(ol[[1]]$beta),
    type = 'l', lwd = lwd, col = col1,
    xlab = xlab, ylab = ylab)
  if (FisherInfo && ('hess' %in% names(ol[[1]]))) {
    # Optionally add dashed lines for Fisher Information "error bars"
    stdb <- 1 / sqrt(diag(ol[[1]]$hess)) 
    lines(x = xvec, y = yfunc(ol[[1]]$beta - stdb), col = col1, lty = 3, lwd = lwd - 1)
    lines(x = xvec, y = yfunc(ol[[1]]$beta + stdb), col = col1, lty = 3, lwd = lwd - 1)
  }
  if (nsols > 1) {
    for (ii in 2:nsols) {
      lines(x = xvec, y = yfunc(ol[[ii]]$beta), lwd = lwd, col = coll[ii])
      if (FisherInfo && ('hess' %in% names(ol[[ii]]))) {
        stdb <- 1 / sqrt(diag(ol[[ii]]$hess)) 
        lines(x = xvec, y = yfunc(ol[[ii]]$beta - stdb), col = coll[ii], lty = 3, lwd = lwd - 1)
        lines(x = xvec, y = yfunc(ol[[ii]]$beta + stdb), col = coll[ii], lty = 3, lwd = lwd - 1)
      }
    }
  }
  return(nsols)
}

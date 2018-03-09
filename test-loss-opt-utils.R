#
# Test / demonstration of codes for estimating smooth univariate curve
# logistic regression (mainly for hazards) with optional Fisher Information
# "uncertainty" indication and enforced "prior base hazard" (also optional) 
#
# Example: Rscript --vanilla test-loss-opt-utils.R
#

source('loss-opt-utils.R')

# generate a lambda(x) test problem to be solved by xentlambda objective
# need to write up foreach enabled CV

lx <- function(x) {
  (1 + sin(x)) * exp(-1 * x * x / pi ^ 2)
}

#
# Model: exponential parameter lx(x); used to generate samples like this:
#
# Prob[Y = 1 | x, l] = 1 - exp(-lx(x) * l) = (def) = yhat
#
# Both the actual probabilities yhat(x), and the outcomes y(x) in {0, 1} are 
# stored; both having a weight l; the data is (y, l, x) or (yhat, l, x) and
# the task is to estimate the function lx(x) from the data.
#

stdx <- pi / 2
lmin <- 0.2
lmax <- 5

n <- 8000
yh.y.x.l <- array(NA, c(n, 4))
for (i in 1:n) {
  l <- runif(1) * (lmax - lmin) + lmin  # random interval / weight
  x <- rnorm(1) * stdx
  yh <- 1 - exp(-lx(x) * l)
  stopifnot(yh >= 0 && yh <= 1)
  y <- ifelse(runif(1) < yh, 1, 0)
  yh.y.x.l[i, ] <- c(yh, y, x, l)
}

ny <- sum(yh.y.x.l[, 2])

# Show what to expect
print(head(yh.y.x.l))
print(tail(yh.y.x.l))
print(sprintf('%i out of %i samples have label = 1', ny, n))

x.range <- c(-3*stdx, 3*stdx)
x.nn <- 60

ell.vec <- 2 ^ seq(from = -7, to = 7, len = 60) #(-7:7)

# Next setup the "inverse problem" to get function lx(x) from the dataset 
# NOTE: with few observations; better results are obtained with true prob labels
# (which are typically unknown)
rep.y <- opt.loss.ell2.pwco.1d(
  yh.y.x.l[, 2],  # col 2: outcome target {0, 1}, col 1: true prob values in [0, 1]
  yh.y.x.l[, 3],
  w = yh.y.x.l[, 4],
  nn = x.nn,
  xrange = x.range,
  ell2.vec = ell.vec, 
  regtyp = 2,
  yreg = -1,  # not in interval [0, 1] and not null means auto-tune (regularize towards baseline)
  lghfunc = lghfunc.xentlambda,
  warm.restart = TRUE)

# Try to select the 'best' parameter lambda.best
rep.y.sel <- auto.pwco.1d.post(rep.y)
lambda.best <- ell.vec[rep.y.sel$idx.min]
print(sprintf('lambda(*) = %e', lambda.best))

#print(rep.y.sel)

rep.y.1 <- opt.loss.ell2.pwco.1d(
  yh.y.x.l[, 2],  # col 2: outcome target {0, 1}, col 1: true prob values in [0, 1]
  yh.y.x.l[, 3],
  w = yh.y.x.l[, 4],
  nn = x.nn,
  xrange = x.range,
  ell2.vec = lambda.best, 
  regtyp = 2,
  yreg = -1,  # not in interval [0, 1] and not null means auto-tune (regularize towards baseline)
  lghfunc = lghfunc.xentlambda,
  warm.restart = TRUE)

# Produce a plot of the solution family (different reg. parameters)

#plot.pwco.1d(rep.y, xrange = x.range, ytyp = 1, dt = 1, FisherInfo = TRUE)

plot.pwco.1d(rep.y.1, xrange = x.range, ytyp = 1, dt = 1, FisherInfo = TRUE)
# next add to plot the true hazard in red
xg <- seq(from = x.range[1], to = x.range[2], len = 200)
lines(x = xg, y = lx(xg), col = 'red', lwd = 4)

# 'rigorous' baseline calculation (using only a single 'bias feature' made of ones)
opt.1 <- opt.loss.ell2(
  lghfunc.xentlambda,
  yh.y.x.l[, 2], 
  as.matrix(array(1, n)), 
  w = yh.y.x.l[, 4],
  ell2 = 0)  # as.matrix(0)

bhzrd.1 <- log(1 + exp(opt.1$beta)) / 1
# 'approximate quick' baseline calculation
bhzrd.2 <- - (1 / mean(yh.y.x.l[, 4])) * log(1 - mean(yh.y.x.l[, 2]))
# and show base hazard as horizontal line
lines(x = c(min(xg), max(xg)), y = rep(bhzrd.1, 2), col = 'black', lwd = 2, lty = 3)
lines(x = c(min(xg), max(xg)), y = rep(bhzrd.2, 2), col = 'black', lwd = 2, lty = 4)

print(c(bhzrd.1, bhzrd.2))

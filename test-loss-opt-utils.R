#
# Test / demonstration of codes for estiamting smooth univariate curve
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

# Show what to expect
print(head(yh.y.x.l))
print(tail(yh.y.x.l))

print(n)
print(sum(yh.y.x.l[, 2]))

x.range <- c(-3*stdx, 3*stdx)
x.nn <- 60

# Next setup the "inverse problem" to get function lx(x) from the dataset 
# NOTE: with few observations; better results are obtained with true prob labels
# (which are typically unknown)
rep.y <- opt.loss.ell2.pwco.1d(
  yh.y.x.l[, 2],  # col 2: outcome target {0, 1}, col 1: true prob values in [0, 1]
  yh.y.x.l[, 3],
  w = yh.y.x.l[, 4],
  nn = x.nn,
  xrange = x.range,
  ell2.vec = c(1/128, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32, 64, 128),
  regtyp = 2,
  yreg = -1,  # not in interval [0, 1] and not null means auto-tune (regularize towards baseline)
  lghfunc = lghfunc.xentlambda,
  warm.restart = TRUE)

# Produce a plot of the solution family (different reg. parameters)
plot.pwco.1d(rep.y, xrange = x.range, ytyp = 1, dt = 1, FisherInfo = TRUE)
# next add to plot the true hazard in red
xg <- seq(from = x.range[1], to = x.range[2], len = 200)
lines(x = xg, y = lx(xg), col = 'red', lwd = 4)
# and show base hazard as horizontal line
bhzrd <- sum(yh.y.x.l[, 2]) / sum(yh.y.x.l[, 4])
lines(x = c(min(xg), max(xg)), y = c(1, 1) * bhzrd, col = 'black', lwd = 2, lty = 3)

# TODO: 
# Next automatically obtain a good ell2 regularization value 
# foreach-based CV example !
#

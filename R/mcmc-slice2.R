## MCMC updates via slice sampling.

## The algorithm is described in detail in
##   Neal R.M. 2003. Slice sampling. Annals of Statistics 31:705-767.
## which describes *why* this algorithm works.  The approach differs
## from normal Metropolis-Hastings algorithms, and from Gibbs
## samplers, but shares the Gibbs sampler property of every update
## being accepted.  Nothing is required except the ability to evaluate
## the function at every point in continuous parameter space.

## Let x0 be the current (possibly multivariate) position in
## continuous parameter space, and let y0 be the probability at that
## point.  To update from (x0, y0) -> (x1, y1), we update each of the
## parameters in turn.  For each parameter
##
##   1. Draw a random number 'z' on Uniform(0, y0) -- the new point
##      must have at least this probability.
##
##   2. Find a region (x.l, x.r) that contains x0[i], such that x.l
##      and x.r are both smaller than z.
##
##   3. Randomly draw a new position from (x.l, x.r).  If this
##      position is greater than z, this is our new position.
##      Otherwise it becomes a new boundary and we repeat this step
##      (the point x1 becomes x.l if x.l < x0[i], and x.r otherwise so
##      that x0[i] is always contained within the interval).

## Because it is generally more convenient to work with log
## probabilities, step 1 is modified so that we draw 'z' by taking
##   y0 - rexp(1)
## All other steps remain unmodified.

## Here, w, lower and upper are vectors
sampler.slice2 <- function(lik, x.init, y.init, w, lower, upper, control) {
  for ( i in seq_along(x.init) ) {
    xy <- slice.1d2(make.unipar(lik, x.init, i),
                   x.init[i], y.init, w[i], lower[i], upper[i])
    x.init[i] <- xy[[1]]
    y.init    <- xy[[2]]
  }

  list(x.init=x.init, y.init=y.init)
}

## Here, w, lower and upper are scalars
slice.1d2 <- function(f, x.init, y.init, w, lower, upper) {
  z <- y.init$loglik - rexp(1)
  r <- slice.isolate2(f, x.init, y.init, z, w, lower, upper)
  slice.sample2(f, x.init, z, r)
}

slice.isolate2 <- function(f, x.init, y.init, z, w, lower, upper) {
  u <- runif(1) * w
  L <- x.init - u
  R <- x.init + (w-u)

  while ( L > lower && f(L)$loglik > z )
    L <- L - w
  while ( R < upper && f(R)$loglik > z )
    R <- R + w

  c(max(L, lower), min(R, upper))
}

slice.sample2 <- function(f, x.init, z, r) {
  r0 <- r[1]
  r1 <- r[2]

  repeat {
    xs <- runif(1, r0, r1)
    ys <- f(xs)
    if ( ys$loglik > z )
      break
    if ( xs < x.init )
      r0 <- xs
    else
      r1 <- xs
  }
  list(xs=xs, ys=ys)
}

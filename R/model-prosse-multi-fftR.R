make.branches.prosse.multi.fftR <- function(control) {
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  dt.max <- control$dt.max
  tc <- control$tc
  h <- control$h

  f.hi <- make.pde.prosse.multi.fftR(nx*r, dx/r, dt.max, 2*h)
  f.lo <- make.pde.prosse.multi.fftR(nx,   dx,   dt.max, 2*h)
  combine.branches.prosse.multi(f.hi, f.lo, control)
}

make.pde.prosse.multi.fftR <- function(nx, dx, dt.max, nd) {
  function(y, len, pars, t0, star) {
    padding <- pars$padding
    ndat <- length(pars$lambda)
    
    nt <- as.integer(ceiling(len / dt.max))
    dt <- len / nt
    b <- pars$b
    mu <- pars$mu
    z <- exp(len * (b-mu))
    Et <- (mu + z*(y[1] - 1)*mu - b*y[1]) / (mu + z*(y[1] - 1)*b - b*y[1])
    
    if ( !(length(y[-1]) %in% (nd * nx)) )
      stop("Wrong size y")
    if ( length(pars$lambda) > (nx-3) )
      stop("Incorrect length pars")
    if ( pars$diffusion <= 0 )
      stop("Invalid diffusion parameter")

    ans <- prosse.integrate.multi.fftR(matrix(y[-1], nx, nd), y[1], star, pars$b, pars$mu, pars$lambda, pars$Qr, pars$Q, pars$drift,
                                 pars$diffusion, nt, dt, nx, ndat, dx, nd,
                                 padding[1], padding[2])
    ## Do the log compensation here, to make the careful calcuations
    ## easier later.
    if ( any(ans < -1e-8) )
      stop("Actual negative D value detected -- calculation failure")
    ans[ans < 0] <- 0
    
    q <- sum(ans) * dx
    ans <- ans / q
    list(log(q),Et, ans)
  }
}

prosse.integrate.multi.fftR <- function(vars, Et, star, b, mu, lambda, Qr, Q, drift, diffusion,
                                  nstep, dt, nx, ndat, dx, nd, nkl, nkr) {
  kern <- fftR.make.kern(-dt * drift, sqrt(dt * diffusion),
                         nx, dx, nkl, nkr)
  fy <- fft(kern)  
  for ( i in seq_len(nstep) ) {
    vars <- fftR.propagate.t.prosse.multi(vars, Et, ii=i, star, b, mu, lambda, Qr, Q, dt, ndat, nd)
    vars <- fftR.propagate.x(vars, nx, fy, nkl, nkr)
  }
  vars
}

fftR.propagate.t.prosse.multi <- function(vars, Et, ii, star, b, mu, lambda, Qr, Q, dt, ndat, nd) {
  i <- seq_len(ndat)
  r <- b - mu
  z <- exp((ii-1) * dt * r)
  tmp1 = mu - b * Et
  tmp2 = z * (Et - 1)
  h <- nd/2
  Et <- (tmp1 + tmp2 * mu) / (tmp1 + tmp2 * b)
  z <- exp(dt * r)
  A <- (z * r * r)/(z * b - mu + (1-z) * b * Et)^2
  if (star==1) {
  	d0 <- A * vars[i,1:h] *exp(-lambda * dt)
    vars[i,1:h] <- t(mapply(function (z) t(Q %*% d0[z,]), i))
    d0 <- A * vars[i,(h+1):(2*h)] *exp(-lambda * dt)
    vars[i,(h+1):(2*h)] <- t(mapply(function (z) t(Q %*% d0[z,]), i))
  } else if (star==2) {
  	dr0 <- A * vars[i,(h+1):(2*h)]
    vars[i,(h+1):(2*h)] <- t(mapply(function (z) t(Qr %*% dr0[z,]), i))
    d0 <- A * vars[i,1:h] * exp(-lambda * dt)
    d0 <- d0 - dr0 * exp(-lambda * dt)
    vars[i,1:h] <- vars[i,(h+1):(2*h)] + t(mapply(function (z) t(Q %*% d0[z,]), i))
  }
  vars
}

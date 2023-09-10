FFTW.ESTIMATE <- -1
FFTW.MEASURE <- 0
FFTW.PATIENT <- 1
FFTW.EXHAUSTIVE <- 2

## 8: branches
## TODO: I am not sure that this will pick up on the edge case of tc
## being *exactly* t0 or t0 + len (and what happens if rounding error
## puts this around a break).
## This might be something to check during the checking phase that tc
## does not come within 1e-8 of a node?
make.branches.prosse.multi.fftC <- function(control) {
  check.fftC()
  
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  dt.max <- control$dt.max
  tc <- control$tc
  flags <- control$flags
  h <- control$h
  
  f.hi <- make.pde.prosse.multi.fftC(nx*r, dx/r, dt.max, 2*h, flags)
  f.lo <- make.pde.prosse.multi.fftC(nx,   dx,   dt.max, 2*h, flags)
  combine.branches.prosse.multi(f.hi, f.lo, control)
}

make.pde.prosse.multi.fftC <- function(nx, dx, dt.max, nd, flags) {
  ptr <- .Call(r_make_prosse_multi_fft, as.integer(nx), as.numeric(dx),
               as.integer(nd), as.integer(flags))
  function(y, len, pars, t0, star, dt=dt.max) {
    nt <- as.integer(ceiling(len / dt.max))
    dt <- len / nt
    b <- pars$b
    mu <- pars$mu
    z <- exp(len * (b-mu))
    Et <- (mu + z*(y[1] - 1)*mu - b*y[1]) / (mu + z*(y[1] - 1)*b - b*y[1])
    
    if ( !(length(y[-1]) %in% (nd * nx)) )
      stop("Wrong size y")
    if ( any(is.na(y)) )
      stop("Illegal NA values in input")
    if ( length(pars$lambda) > (nx-3) )
      stop("Incorrect length pars")
    if ( pars$diffusion <= 0 )
      stop("Invalid diffusion parameter")
    
    ans <- .Call(r_do_integrate_prosse_multi,
                 ptr, y[-1], y[1], as.integer(star), pars$b, pars$mu, pars$lambda,
                 as.numeric(pars$Q), as.numeric(pars$Qr), pars$drift, pars$diffusion,
                 nt, dt, pars$padding)

    ## Do the log compensation here, to make the careful calcuations
    ## easier later.
    if ( any(ans < -1e-8) )
      stop("Actual negative D value detected -- calculation failure")
    ans[ans < 0] <- 0
    
    q <- sum(ans) * dx
    ans <- ans / q
    list(log(q), Et, ans)
  }
}


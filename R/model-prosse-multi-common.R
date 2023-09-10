exp.x <- function(x,r) {
	r*exp(x)
}

starting.point.prosse.multi <- function(tree, lik, q.div=5, states, states.sd=NULL) {
  nq <- length(get.cache(lik)$args$q)
  pars.bd <- starting.point.bd(tree)
  if ( pars.bd[1] > pars.bd[2] )
    p <- c(pars.bd,pars.bd[1],rep(pars.bd[1]/q.div,nq), rep(0.5,nq)) #b, mu, lambda, q, p
  else
    p <- c(pars.bd[1],pars.bd[1]/q.div,pars.bd[1],rep(pars.bd[1]/q.div,nq), rep(0.5,nq))
    
  lik.bm <- make.bm(tree, states, states.sd,
                    control=list(method="pruning", backend="C"))
  s2 <- max(stats::coef(find.mle(lik.bm, .1)),10^-3)
  c(p, 0, diffusion=s2)
}

## I use a number of elements of pars.
## pars[[i]]${b,mu,lambda,drift,diffusion,padding}
## pars$tr
expand.pars.prosse.multi <- function(lambda, args, ext, pars) {
  pars.use <- vector("list", 2)
  for ( i in c(1,2) ) {
    x <- list()
    pars.use[[i]] <-
      list(x=ext$x[[i]], # May screw other things up (was $x[i])
           lambda=do.call(lambda, c(ext$x[i], pars[args$lambda])),
           drift=pars[args$drift],
           diffusion=pars[args$diffusion],
           padding=ext$padding[i,],
           ndat=ext$ndat[i],
           nx=ext$nx[i])
  }
  names(pars.use) <- c("hi", "lo")
  pars.use$tr <- ext$tr
  pars.use
}

make.pars.prosse.multi <- function(cache) {
  args <- cache$args
  dt.max <- cache$control$dt.max
  k <- cache$control$k
  h <- cache$control$h

  function(pars) {
    names(pars) <- NULL # Because of use of do.call, strip names
	
	b <- pars[args$b]
    mu <- pars[args$mu]
	q <- pars[args$q]
	p <- pars[args$p]
    drift <- pars[args$drift]
    diffusion <- pars[args$diffusion]

    ext <- quasse.extent(cache$control, drift, diffusion)
    ## This would confirm the translation:
    ##   all.equal(ext$x[[1]][ext$tr], ext$x[[2]])

    ## Parameters, expanded onto the extent:
    pars <- expand.pars.prosse.multi(cache$lambda, args, ext, pars)

    check.pars.prosse.multi(b, mu, pars$hi$lambda, q, p, drift, diffusion)
    
    #calculate q_r matrix
    pars$nam <- expand.grid(lapply(k,function (i) 1:i))
    q.r <- diag(dim(pars$nam)[1])
    idx <- cbind(rep((1:h), each=k[1]-1), rep(unlist(lapply(1:k[1], function(i) (1:k[1])[-i])), h/k[1]) + rep(seq(0,h-k[1],k[1]),each=k[1]*(k[1]-1)))
    q.r[idx] <- rep(q[1:(k[1]*(k[1]-1))], h/k[1])
    if (length(k)>1) {
    for (i in 2:length(k)) {
    	a <- h/prod(k[1:(i-1)])
    	idx <- cbind(rep((1:a), each=k[i]-1), rep(unlist(lapply(1:k[i], function(j) (1:k[i])[-j])), a/k[i]) + rep(seq(0,a-k[i],k[i]),each=k[i]*(k[i]-1)))
    	tmp <- idx
    	for (j in 1:(max(idx)-1)) {
    		tmp[idx==(j+1)] <- tmp[idx==(j+1)] + j*(prod(k[1:(i-1)])-1)
    	}
    	a <- prod(k[1:(i-1)])-1
    	idx <- cbind(rep(tmp[,1],a+1)+rep(0:a,each=dim(tmp)[1]),rep(tmp[,2],a+1)+rep(0:a,each=dim(tmp)[1]))
    	q.r[idx] <- rep(q[sum((k*(k-1))[1:(i-1)]) + (1:(k[i]*(k[i]-1)))], h/k[i])
    }
    }
    diag(q.r) <- 0
    diag(q.r) <- -rowSums(q.r)
    
    #calculate q matrix
    q.m <- diag(dim(pars$nam)[1])
    idx <- cbind(rep((1:h), each=k[1]-1), rep(unlist(lapply(1:k[1], function(i) (1:k[1])[-i])), h/k[1]) + rep(seq(0,h-k[1],k[1]),each=k[1]*(k[1]-1)))
    q.m[idx] <- rep(q[1:(k[1]*(k[1]-1))]*(1-p[1:(k[1]*(k[1]-1))]), h/k[1])
    if (length(k)>1) {
    for (i in 2:length(k)) {
    	a <- h/prod(k[1:(i-1)])
    	idx <- cbind(rep((1:a), each=k[i]-1), rep(unlist(lapply(1:k[i], function(j) (1:k[i])[-j])), a/k[i]) + rep(seq(0,a-k[i],k[i]),each=k[i]*(k[i]-1)))
    	tmp <- idx
    	for (j in 1:(max(idx)-1)) {
    		tmp[idx==(j+1)] <- tmp[idx==(j+1)] + j*(prod(k[1:(i-1)])-1)
    	}
    	a <- prod(k[1:(i-1)])-1
    	idx <- cbind(rep(tmp[,1],a+1)+rep(0:a,each=dim(tmp)[1]),rep(tmp[,2],a+1)+rep(0:a,each=dim(tmp)[1]))
    	q.m[idx] <- rep(q[sum((k*(k-1))[1:(i-1)]) + (1:(k[i]*(k[i]-1)))]*(1-p[sum((k*(k-1))[1:(i-1)]) + (1:(k[i]*(k[i]-1)))]), h/k[i])
    }
    }
    diag(q.m) <- diag(q.r)

    #calculate e^(Q*dt)
    pars$lo$Qr <- pars$hi$Qr <- as.matrix(Matrix::expm(q.r*dt.max))
	pars$lo$Q <- pars$hi$Q <- as.matrix(Matrix::expm(q.m*dt.max))
	pars$lo$b <- pars$hi$b <- b
	pars$lo$mu <- pars$hi$mu <- mu
	pars$lo$drift <- pars$hi$drift <- drift
	pars$lo$diffusion <- pars$hi$diffusion <- diffusion

    pars
  }
}

combine.branches.prosse.multi <- function(f.hi, f.lo, control) {
  nx <- control$nx
  dx <- control$dx
  tc <- control$tc
  r <- control$r
  h <- control$h
  eps <- log(control$eps)
  dt.max <- control$dt.max

  ## This is hacky version of the log compensation.  It also reduces
  ## the stepsize when bisecting a branch.  It doesn't seem to change
  ## much on 
  careful <- function(f, y, len, pars, t0, star, dt.max) {
    ans <- f(y, len, pars, t0, star)
    if ( ans[[1]] > eps ) { # OK
      ans
    } else {
      if ( control$method == "fftC" ||
           control$method == "fftR" )
        dt.max <- dt.max / 2 # Possibly needed
      len2 <- len/2
      ans1 <- Recall(f, y,         len2, pars, t0,        star, dt.max)
      ans2 <- Recall(f, c(ans1[[2]],as.numeric(ans1[[3]])), len2, pars, t0 + len2, star, dt.max)
      ans2[[1]][[1]] <- ans1[[1]][[1]] + ans2[[1]][[1]]
      ans2
    }
  }

  ## Start by normalising the input so that eps up there make
  ## sense...
  function(y, len, pars, t0, star, idx) {
    if ( t0 < tc ) {
      dx0 <- dx / r
      nx0 <- nx * r
    } else {
      dx0 <- dx
      nx0 <- nx
    }

    ## Here, we also squash all negative numbers.
    #if ( any(y < -1e-8) )
    #  stop("Actual negative D value detected -- calculation failure")
    #y[y < 0] <- 0
    #q0 <- sum(y[-1]) * dx0
    #if ( q0 <= 0 )
    #  stop("No positive D values")
    #y[-1] <- y[-1] / q0
    #lq0 <- log(q0)
    lq0 <- 0
    
    if ( t0 >= tc ) {
      ans <- careful(f.lo, y, len, pars$lo, t0, star, dt.max)
    } else if ( t0 + len < tc ) {
      ans <- careful(f.hi, y, len, pars$hi, t0, star, dt.max)
    } else {
      len.hi <- tc - t0
      ans.hi <- careful(f.hi, y, len.hi, pars$hi, t0, star, dt.max)

      y.lo <- ans.hi[[3]][pars$tr,]
      lq0 <- lq0 + ans.hi[[1]]
      if ( nrow(y.lo) < nx )
        y.lo <- rbind(y.lo, matrix(0, nx - length(pars$tr), 2*h))
	  
	  y.lo <- c(ans.hi[[2]],as.numeric(y.lo))
      ## Fininshing up with the low resolution branch...
      ans <- careful(f.lo, y.lo, len - len.hi, pars$lo, tc, star, dt.max)
    }

    list(ans[[1]] + lq0, ans[[2]], ans[[3]])
  }
}

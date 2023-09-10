## Models should provide:
##   1. make
##   2. info
##   3. make.cache, including initial tip conditions
##   4. initial.conditions(init, pars,t, idx)
##   5. rootfunc(res, pars, ...)

## Common other functions include:
##   stationary.freq
##   starting.point
##   branches

## 1: make
make.prosse.multi <- function (tree, traits, states, states.sd, lambda, control=NULL) {
	cache <- make.cache.prosse.multi(tree, traits, states, states.sd, lambda, control)
	all.branches <- make.all.branches.prosse.multi(cache, cache$control)
	rootfunc <- make.rootfunc.prosse.multi(cache)
	f.pars <- make.pars.prosse.multi(cache)
    ll <- function (pars, condition.surv=TRUE, root=ROOT.OBS, root.f=NULL) {
		pars2 <- f.pars(pars)
		ans <- all.branches(pars2)
        rootfunc(ans, pars2, condition.surv, root, root.f)
	}
	class(ll) <- c("prosse","dtlik","function")
	ll
}

## 2: info
make.info.prosse.multi <- function(k, lambda, phy) {
  argnames <- default.argnames.prosse.multi(k, lambda)
  list(name="prosse_multi",
       ## Parameters:
       np=NA,
       argnames=argnames,
       ## Variables:
       ny=NA,
       k=NA,
       idx.e=NA,
       idx.d=NA,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=FALSE, # not for many models.
       ## These are optional
       doc=NULL,
       reference=NA,
       ## These are special to ProSSE:
       lambda=lambda
	)
}
default.argnames.prosse.multi <- function(k, lambda) {
  c("b","mu",
    sprintf("l.%s", names(formals(lambda))[-1]),
    unlist(lapply(1:length(k), function (j) sprintf("q%s%s.%s", rep(as.character(1:k[j]), each=k[j]-1), unlist(lapply(1:k[j], function(i) as.character(1:k[j])[-i])),j))),
    unlist(lapply(1:length(k), function (j) sprintf("p%s%s.%s", rep(as.character(1:k[j]), each=k[j]-1), unlist(lapply(1:k[j], function(i) as.character(1:k[j])[-i])),j))),
    "drift","diffusion"
   )
}

## 3: make.cache (& initial conditions)
make.cache.prosse.multi <- function (tree, traits, states, states.sd, lambda, control) {
	## 1: tree
  	#tree <- check.tree(tree)
  	
  	## 2: states & errors
  	tmp <- check.states.quasse(tree, states, states.sd)
  	states <- tmp$states
 	states.sd <- tmp$states.sd
 	
 	## 3: Control structure (lots of checking!)
 	control <- check.control.quasse(control, tree, states)
 	
	cache <- make.cache.tree.prosse(tree)
	cache$states  <- states
    cache$states.sd <- states.sd
	cache$control <- control
	
	## 4: Speciation/extinction functions
    n.lambda <- check.f.quasse(lambda)
    cache$lambda <- lambda
    
    ## 5: trait parameters
    k <- sapply(1:ncol(traits),function (i) length(unique(traits[,i])))
    n.q <- sum(k*(k-1))
    n.p <- n.q
    cache$traits <- as.matrix(traits)
    cache$control$k <- k
    cache$control$h <- prod(k)
    
    n.args   <- 2 + n.lambda + n.q + n.p + 2
    args <- list(b = 1,
                 mu = 2,
                 lambda = 2 + seq_len(n.lambda),
                 q = 2 + n.lambda + seq_len(n.q),
                 p = 2 + n.lambda + n.q + seq_len(n.p),
                 drift = 2 + n.lambda + n.q + n.p + 1,
                 diffusion = 2 + n.lambda + n.q + n.p + 2)
    cache$args <- args
	
	cache$info <- make.info.prosse.multi(k, lambda, tree)
	
	cache
}

## 4: initial.conditions
initial.tip.prosse.multi <- function (cache, control, x) {
	nx <- control$nx * control$r
	h <- control$h
    k <- control$k
    npad <- nx - length(x)
    states <- cache$states
    states.sd <- cache$states.sd
    traits <- cache$traits
	group <- cache$group # a list with each cell corresponding to a clade, which records the order of branches to calculate within the clade
	n.tip <- cache$n.tip # total number of tips in the tree
	len <- cache$len # branch length of each branch
	#init is a function that generates y=all possible scenarios of the tip initial states in a clade, target=idx of tip branches, t=tip branch length
	init <- function (group2, n.tip, len, nx, npad, states, states.sd, traits, k, h) {
		target <- group2[group2<=n.tip] # find the tip branches
		# the initial state is y, which is a matrix with rows are tips, columns are c(E,DI,DR)
		if (length(target)==1) {
            if (length(k)>1) {
                trait.idx <- sum(sapply(1:(length(k)-1), function (i) (traits[target,i]-1)*prod(k[-c(1:i)])))+traits[target,length(k)]
            } else {
                trait.idx <- traits[target,1]
            }
			y <- c(0, #E
			       rep(0,nx*h), #DI
				   rep(0,nx*(trait.idx-1)), #DR
				   dnorm(x, states[target], states.sd[target]), rep(0, npad),
				   rep(0,nx*(h-trait.idx)))
            list(target=target,t=len[target],y=y)
        } else if (length(target)>1) {
            if (length(k)>1) {
                trait.idx <- rowSums(sapply(1:(length(k)-1), function (i) (traits[target,i]-1)*prod(k[-c(1:i)])))+traits[target,length(k)]
            } else {
                trait.idx <- traits[target,1]
            }
			y <- mapply(function(mean, sd, t.idx) 
					c(0, #E
				   	rep(0,nx*(t.idx-1)), #DI
				   	dnorm(x, mean, sd), rep(0, npad),
				  	rep(0,nx*(h-t.idx)),
				  	rep(0,nx*h)), #DR
                    states[target], states.sd[target], trait.idx, SIMPLIFY=FALSE)
            y[[1]] <- c(0,
                   rep(0,nx*h),
				   rep(0,nx*(trait.idx[1]-1)),
				   dnorm(x, states[target[1]], states.sd[target[1]]), rep(0, npad),
				   rep(0,nx*(h-trait.idx[1])))        
            list(target=target,t=len[target],y=y)
		} else {
			list(NULL)
		}
	}
	lapply(group,init,n.tip,len, nx, npad, states, states.sd, traits, k, h) #apply init to each cell in the group list.
}

make.initial.conditions.prosse.multi <- function(control) {
  h <- control$h
  tc <- control$tc
  r <- control$r
  nx.lo <- control$nx
  nx.hi <- nx.lo * r
  dx <- control$dx

  ## There is the chance that we could be slightly off on the depth
  ## by rounding error.  Because of this, I've done the testing
  ## against the *length* of the data, and then checked that the time
  ## is appropriate (to within eps of the correct value).  It is
  ## possible that two different branches with different numbers of
  ## nodes that finish right at the critical interval might have
  ## conflicting lengths.
  eps <- 1e-8
  function(init, pars, t, idx) {
    if ( length(init[[1]]) != length(init[[2]]) )
      stop("Data have incompatible length")
    b <- pars$lo$b
    
    if ( t < tc ) {
        nx <- nx.hi
        lambda <- pars$hi$lambda
        dx <- dx/r
    } else {
        nx <- nx.lo
        lambda <- pars$lo$lambda
    }
    
    ndat <- length(lambda)
    i <- seq_len(ndat)
    
	#DI
    res <- init[[1]][i,1:h] * init[[2]][i,1:h] * b * dx
    
    #DR
    res <- cbind(res,(init[[1]][i,(h+1):(2*h)] * init[[2]][i,1:h] + init[[2]][i,(h+1):(2*h)] * init[[1]][i,1:h]) * b * dx)
    res <- rbind(res, matrix(as.integer(0), nx - ndat, 2*h))
    
    q <- sum(res) * dx
    res <- res / q
    
    c(log(q),as.numeric(res))
  }
}

make.rootfunc.prosse.multi <- function(cache) {
  root.idx <- cache$root
  nx <- cache$control$nx
  dx <- cache$control$dx
  h <- cache$control$h
  root.t <- max(cache$depth)
  
  function(res, pars, condition.surv, root, root.f, intermediates) {
    lq <- res$lq
    e.root <- res$e.root
    d.root <- res$vals[seq_len(pars$lo$ndat),(h+1):(2*h)]
	
    func <- function (t,y,pars) {
        b <- pars[1]
        mu <- pars[2]
        lambda <- pars[3]
        Ei <- y[1]
        Er <- y[2]
        list(c(-(b+mu)*Ei + b*Ei*Ei + mu,
        -(b+mu+lambda)*Er + b*Er*Er + mu + lambda * Ei))
    }

    er.root <- sapply(pars$lo$lambda, function (i) ode(y=c(0,1),times=c(0,root.t),func=func,parms=c(pars$lo$b,pars$lo$mu,i))[2,2])

    root.p.prosse <- function (d.root, pars, e.root, er.root, root, root.f, condition.surv) {
		if ( condition.surv ) {
      		b <- pars$lo$b
            d.root <- d.root / (b * (1 - e.root) * (1- er.root))
    	}
        root.p <- root.p.quasse(d.root, pars$lo, root, root.f)
    	sum(root.p * d.root) * dx
	}
    d.root <- apply(d.root, 2, root.p.prosse, pars, e.root, er.root, root, root.f, condition.surv)
    root.p <- root.p.calc(d.root, pars, root)
	
	loglik <- log(sum(root.p * d.root)) + sum(lq)

    loglik
  }
}

######################################################################
## Extra core stuff:
#this is similar to make.all.branches.dtlik
make.all.branches.prosse.multi <- function (cache,control) {
	branches <- make.branches.prosse.multi(cache,control)
	initial.conditions <- make.initial.conditions.prosse.multi(control)
	function (pars) {
		cache$y <- initial.tip.prosse.multi(cache, control, pars[[1]]$x)
		all.branches.matrix.prosse.multi(pars,cache,initial.conditions,branches)
	}
}

all.branches.matrix.prosse.multi <- function (pars,cache,initial.conditions,branches) {
	len <- cache$len #branch length
  	depth <- cache$depth #branch starting time
  	children <- cache$children #children of each node
  	order <- cache$order #parent of each node
  	root <- cache$root #root label
	group <- cache$group #list of clades
	stars <- cache$stars #indicator for which ODE to use
	h <- cache$control$h
	
  	n <- length(len)
  	lq <- rep(0, n)  #this is the log of DI+DR after the calculation on each branch. It is used to standardize DI and DR to avoid rounding error due to extremely small likelihood value
 	Et <- rep(0,n)
 	n.tip <- cache$n.tip

  	y <- cache$y #initial tip states
  	branch.init <- branch.base <- vector("list", n)
	
	#start from each clade
    for ( x in 1:length(y) ) {
        idx <- y[[x]]$target #tip index in each clade
        if (length(idx)==1) { #for clade with single tip
        	branch.init[idx] <- list(y[[x]]$y)
        	ans <- branches(y[[x]]$y, y[[x]]$t, pars, 0, star=2)
        	lq[idx] <- ans[[1]]
        	Et[idx] <- ans[[2]]
        	branch.base[[idx]] <- ans[[3]]
        } else if (length(idx)>1){ #for clade more than one tip
        	idx2 <- group[[x]]
        	idx2 <- idx2[idx2>n.tip]
        	#calculate along tip branch
        	for (j in 1:length(idx)) { 
        		branch.init[idx[j]] <- y[[x]]$y[j]
        		ans <- branches(y[[x]]$y[[j]], y[[x]]$t[j], pars, 0, star=1) #use the ODE for DI*, so star=1
         		lq[idx[j]] <- ans[[1]]
         		Et[idx[j]] <- ans[[2]]
        		branch.base[[idx[j]]] <- ans[[3]]
        	}
        	#calculate along internal branch
         	for (j in idx2[-length(idx2)]) {
         		tmp <- children[j,!is.element(children[j,],group[[x]])] 
         		if (length(tmp)>0) {
         			branch.base[[tmp]][,(h+1):(2*h)] <- 0
         		} #find children that leads to another species and set its DR=0
        		y.in <- initial.conditions(branch.base[children[j,]], pars, depth[j]) #calculate the boundary condition of the branch
        		if ( !all(is.finite(y.in)) )
      				stop("Bad initial conditions: calculation failure along branches?")
    			branch.init[[j]] <- y.in[-1]
    			ans <- branches(c(Et[children[j,1]], y.in[-1]), len[j], pars, depth[j], star=1)
                lq[j] <- y.in[1] + ans[[1]]
    			Et[j] <- ans[[2]]
    			branch.base[[j]] <- ans[[3]]
        	}
        	#calculate the boundary conditions for the basal branch
        	j <- idx2[length(idx2)]
        	tmp <- children[j,!is.element(children[j,],group[[x]])] 
         	if (length(tmp)>0) {
         		branch.base[[tmp]][,(h+1):(2*h)] <- 0
         	} #find children that leads to another species and set its DR=0
        	y.in <- initial.conditions(branch.base[children[j,]], pars, depth[j])
        	if ( !all(is.finite(y.in)) )
      			stop("Bad initial conditions: calculation failure along branches?")
    		branch.init[[j]] <- y.in[-1]
            lq[j] <- y.in[1]
    		#calculate along the basal branch
    		if (j!=root) {
        		ans <- branches(c(Et[children[j,1]], y.in[-1]), len[j], pars, depth[j], star=2) #use the ODE for DI, so star=2
        		lq[j] <- lq[j] + ans[[1]]
    			Et[j] <- ans[[2]]
    			branch.base[[j]] <- ans[[3]]
        	}
        } else { #if the group is not a clade, but the branch linking a clade to another clade that it is nested in, then we need to calculate this branch first.
        	for (j in group[[x]]) {
        		y.in <- initial.conditions(branch.base[children[j,]], pars, depth[j])
        		if ( !all(is.finite(y.in)) )
      				stop("Bad initial conditions: calculation failure along branches?")
    			branch.init[[j]] <- y.in[-1]
                lq[j] <- y.in[1]
    			if (j!=root) {
    				ans <- branches(c(Et[children[j,1]], y.in[-1]), len[j], pars, depth[j], star=2)
    				lq[j] <- lq[j] + ans[[1]]
    				Et[j] <- ans[[2]]
    				branch.base[[j]] <- ans[[3]]
        		}
        	}
        }
    }

    for ( i in order ) {
    	if (i!=root) {
    	y.in <- initial.conditions(branch.base[children[i,]], pars, depth[i])
    	if ( !all(is.finite(y.in)) )
      		stop("Bad initial conditions: calculation failure along branches?")
    	branch.init[[i]] <- y.in[-1]
    	ans <- branches(c(Et[children[i,1]], y.in[-1]), len[i], pars, depth[i], star=2)
        lq[i] <- y.in[1] + ans[[1]]
    	Et[i] <- ans[[2]]
    	branch.base[[i]] <- ans[[3]]
    	}
  	}
 	
  	y.in <- initial.conditions(branch.base[children[root,]], pars, depth[root])
  	branch.init[[root]] <- y.in[-1]
  	ans <- branches(c(Et[children[root,1]], y.in[-1]), len[root], pars, depth[root], star=2)
    lq[root] <- y.in[1] + ans[[1]]
    Et[root] <- ans[[2]]
    branch.base[[root]] <- ans[[3]]
  	
  list(init=branch.init, base=branch.base, lq=lq, vals=branch.base[[root]], e.root=Et[root])
}

#this is similar to make.branches.dtlik, except we have two suits of ODEs, which ODE to use is determined by star, with 1=DI* and DR, 2=DI and DR.
make.branches.prosse.multi <- function(cache,control) {
	if ( control$method == "fftC" ) {
    	branches <- make.branches.prosse.multi.fftC(control)
  	} else if ( control$method == "fftR" ) {
    	branches <- make.branches.prosse.multi.fftR(control)
	}
}
	
check.pars.prosse.multi <- function(b, mu, lambda.x, q, p, drift, diffusion) {
  if ( any(!is.finite(c(b, mu, lambda.x, q, p, drift, diffusion))) )
    stop("Non-finite/NA parameters")
  if ( b < 0 || mu < 0 || any(lambda.x < 0) || any(q < 0) || any(p < 0) || diffusion <= 0 )
    stop("Illegal negative parameters")
  if (any(p > 1))
  	stop("p cannot exceed 1")
  if ( b == 0 )
    stop("No positive b; cannot compute likelihood")
}

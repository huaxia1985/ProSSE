make.prosse.lambda.state <- function (tree, types, states, states.sd, lambda, control=NULL) {
	cache <- make.cache.prosse.lambda.state(tree, types, states, states.sd, lambda, control)
	f.pars <- make.pars.prosse.lambda.state(cache)
	all.branches <- make.all.branches.prosse.lambda.state(cache)
	rootfunc <- function (res, pars, t=max(cache$depth), condition.surv) {
  		vals <- res$vals
  		lq <- res$lq
  		d.root <- vals[3] #the root must be in state R
  
  		if ( condition.surv ) { #condition the likelihood on that the survivial of the root
    		b <- pars[1]
    		e.root <- vals[1]
    		func <- function (t,y,pars) {
    			b <- pars[1]
    			mu <- pars[2]
    			lambda <- pars[3]
    			Ei <- y[1]
    			Er <- y[2]
    			list(c(-(b+mu)*Ei + b*Ei*Ei + mu,
    			-(b+mu+lambda)*Er + b*Er*Er + mu + lambda * Ei))
    		}
    		er.root <- ode(y=c(0,1),times=c(0,t),func=func,pars)[2,2]
    		d.root <- d.root / (b * (1 - e.root) * (1 - er.root))
  		}

  		loglik <- log(d.root) + sum(lq)

  		loglik
  	}
	ll <- function (pars, condition.surv=T) {
		pars2 <- f.pars(pars)
		ans <- all.branches(pars2)
		rootfunc(ans, pars, max(cache$depth),condition.surv)
	}
	class(ll) <- c("prosse","dtlik","function")
	ll
}

## 2: info
make.info.prosse.lambda.state <- function (lambda, phy) {
	list(name="prosse_lambda_state",
		 ## Parameters:
		 np=NA,
		 argnames=default.argnames.prosse.lambda.state(lambda),
		 ## Variables:
		 ny=NA,
		 k=NA,
		 idx.e=NA,
		 idx.d=NA,
		 ## Phylogeny:
		 phy=phy,
		 ## Inference:
		 ml.default="subplex",
		 mcmc.lowerzero=TRUE)
	}

default.argnames.prosse.lambda.state <- function(lambda) {
  c("b","mu","lambda.m",
    sprintf("l.%s", names(formals(lambda))[-1])
  )
}

## 3: make.cache (& initial conditions)
make.cache.prosse.lambda.state <- function (tree, types, states, states.sd, lambda, control) {
  	#tree <- check.tree(tree)
	cache <- make.cache.tree.prosse.lambda(tree,types)
    #parameters
    n.lambda <- diversitree:::check.f.quasse(lambda)
    n.args   <- 3 + n.lambda
    args <- list(b = 1,
                 mu = 2,
                 lambda.m = 3,
                 lambda.c = 3 + seq_len(n.lambda)
                 )
    cache$args <- args
    #lambda
    control <- check.control.prosse.lambda.state(control, states)
    cache$nx <- nx <- control$nx
  	cache$dx <- dx <- control$dx
  	xmid <- control$xmid
  	cache$x <- seq(xmid - dx*ceiling((nx - 1)/2), length.out=nx, by=dx)
  	type <- seq_len(max(types))
    cache$states <- sapply(type,function (i) states[which(types==i)[1]])
    cache$states.sd <- sapply(type,function (i) states.sd[which(types==i)[1]])
	#info
   	cache$info <- make.info.prosse.lambda.state(lambda, tree)
   	cache$h <- 1
	cache
}

check.control.prosse.lambda.state <- function(control, states) {
  xr <- diff(range(states))
  xr.mult <- if ( "xr.mult" %in% names(control) )
    control$xr.mult else 5
  defaults <- list(nx=1024,
                   dx=xr * xr.mult / 1024,
                   xmid=mean(range(states))
                   )
  nx.changed <- "nx" %in% names(control)
  dx.changed <- "dx" %in% names(control)
  control <- if ( is.null(control) ) {
    	defaults
    } else {
    	modifyList(defaults, control)
    }
  if ( dx.changed && !nx.changed ) {
    control$nx <- 2^ceiling(log2(xr * xr.mult / control$dx))
  } else if ( nx.changed && !dx.changed ) {
    control$dx <- xr * xr.mult / control$nx
  }
  
  control
}

## 4: make.pars
make.pars.prosse.lambda.state <- function (cache) {
	args <- cache$args
	x <- cache$x
	function (pars) {
		b <- pars[args$b]
   		mu <- pars[args$mu]
   		lambda.m <- pars[args$lambda.m]
		lambda.c <- lambda(x, pars[args$lambda.c])
		names(pars) <- NULL
		list(b=b,mu=mu,lambda.m=lambda.m,lambda.c=lambda.c)
	}
}

#5. initial conditions
initial.tip.prosse.lambda.state <- function (cache) {
    nx <- cache$nx
    types <- cache$types
    traits <- cache$traits
	group <- cache$group # a list with each cell corresponding to a clade, which records the order of branches to calculate within the clade
	n.tip <- cache$n.tip # total number of tips in the tree
	len <- cache$len # branch length of each branch
	#init is a function that generates y=all possible scenarios of the tip initial states in a clade, target=idx of tip branches, t=tip branch length
	init <- function (group2,n.tip,len) {
		target <- group2[group2<=n.tip] # find the tip branches
		# the initial state is y, which is a matrix with rows are tips, columns are c(E,DI,DR)
		k <- length(target)
		if (k==1) {
            type <- types[target]
            if (type > 0) {
				y <- c(0, #E
			       rep(0,nx), #DI
				   rep(1,nx)  #DR
				   )
			} else {
				y <- c(0, #E
			       0, #DI
				   1  #DR
				   )
			}
            list(target=target,t=len[target],y=y)
        } else if (k>1) {
            type <- types[target[1]]
            if (type > 0) {
				y <- rep(list(
						c(0, #E
				   		rep(1,nx),  #DI
				  		rep(0,nx))  #DR
                    	),k)
           		y[[1]] <- c(0,
                   rep(0,nx),
				   rep(1,nx))
			} else {
				y <- rep(list(c(0,1,0)),k)
				y[[1]] <- c(0,0,1)
			}		          
            list(target=target,t=len[target],y=y)
		} else {
			list(NULL)
		}
	}
	lapply(group,init,n.tip,len) #apply init to each cell in the group list.
}

make.all.branches.prosse.lambda.state <- function (cache) {
	nx <- cache$nx

	branches <- function(y, len, pars, t0, star, type) {
		i <- seq_len(nx)
		
		b <- pars$b
		mu <- pars$mu
		r <- b-mu
		z <- exp(len * r)
		e0 <- y[1]
  		y[1] <- (mu + z*(e0 - 1)*mu - b*e0) / (mu + z*(e0 - 1)*b - b*e0)
  		A <- (z * r * r)/(z * b - mu + (1-z)*b*e0)^2
  		if (type==0) {
  			lambda <- pars$lambda.m
  			z2 <- exp(- len * lambda)
  			if (star==1) {
  				y[-1] <- A * y[-1] * z2
  			} else if (star==2) {
  				dr0 <- y[3]
  				y[3] <- A * y[3]
  				y[2] <- y[3] + (y[2]-dr0) * A * z2
    		}
  		} else {
  			lambda <- pars$lambda.c
  			z2 <- exp(- len * lambda)
  			vars <- matrix(y[-1],nx,2)
  			if (star==1) {
    			vars[i,1] <- A * vars[i,1] * z2
    			vars[i,2] <- A * vars[i,2] * z2
    		} else if (star==2) {
    			dr0 <- vars[i,2]
    			vars[i,2] <- A * vars[i,2]
    			vars[i,1] <- vars[i,2] + (vars[i,1]-dr0) * A * z2
    		}
    		y <- c(y[1],as.numeric(vars))
  		}
  		if (any((y[-1] <= 10^-5) * (y[-1]>0)==1)) {
  			lq <- sum(y[-1])
  			y[-1] <- y[-1] / lq
  			lq <- log(lq)
  		} else {
  			lq <- 0
  		}
  		c(lq, y) 	
		}
	function (pars) {
		cache$y <- initial.tip.prosse.lambda.state(cache)
		all.branches.matrix.prosse.trait.lambda.state(pars,cache,initial.conditions=initial.conditions.prosse.trait.lambda.state,branches)
	}
}


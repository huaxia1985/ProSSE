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
make.prosse.sp <- function (tree,control=list(compiled=FALSE,backend="gslode",condition.surv=T)) {
	cache <- make.cache.prosse.sp(tree,control)
	all.branches <- make.all.branches.prosse.sp(cache,cache$control)
	rootfunc <- function (res, pars, condition.surv) {
  		vals <- res$vals
  		lq <- res$lq
  		d.root <- vals[3:4] #the root must be in state R
  		er.root <- res$er.root
  		
  		if ( condition.surv ) { #condition the likelihood on that the survivial of the root
    		b <- pars[1]
    		e.root <- vals[1]
    		d.root[1] <- d.root[1] / (b * (1 - e.root) * (1 - er.root))
    		d.root[2] <- d.root[2] / (b * (1 - er.root) * (1 - er.root))
  		}

  		loglik <- log(sum(d.root)) + sum(lq)

  		loglik
  	}
	ll <- function (pars,condition.surv=control$condition.surv) {
		ans <- all.branches(pars)
		rootfunc(ans, pars, condition.surv)
	}
	class(ll) <- c("prosse","dtlik","function")
	ll
}

## 3: make.cache (& initial conditions)
make.cache.prosse.sp <- function (tree,control) {
	## 1: tree
  	#tree <- check.tree(tree)
  	
  	## 2: Control structure
	cache <- make.cache(tree)
	cache$len[cache$root] <- tree$root.depth-max(cache$depth)
	control <- check.control.ode(control)
	cache$control <- control
	
	cache$info <- make.info.prosse.sp(tree)
	
	cache
}

initial.conditions.prosse.sp <- function(init, pars, t)
  c(init[[1]][1],init[[1]][2],
    (init[[1]][3] * init[[2]][5] + init[[1]][5] * init[[2]][3]) * pars[1],
    (init[[1]][4] * init[[2]][6] + init[[1]][6] * init[[2]][4]) * pars[1],
    init[[1]][5] * init[[2]][5] * pars[1],
    init[[1]][6] * init[[2]][6] * pars[1])

######################################################################
## Extra core stuff:
#this is similar to make.all.branches.dtlik
make.all.branches.prosse.sp <- function (cache,control) {
	branches <- make.branches.prosse.sp(cache$info,control)
	function (pars) {
		cache$y <- rep(list(c(0,1,1,0,0,0)),cache$n.tip)
		all.branches.matrix.prosse.sp(pars,cache,initial.conditions=initial.conditions.prosse.sp,branches)
	}
}

make.branches.prosse.sp <- function(info,control) {
	eps <- control$eps
	ode <- make.ode(info,control) #make a list of ODEs, with the first cell for DI* and the second cell for DI.
	branches <- function(y, len, pars, t0) {ode(y, c(t0, t0+len), pars)}
  	make.branches.comp.sp(branches, eps)
}

all.branches.matrix.prosse.sp <- function (pars,cache,initial.conditions,branches) {
	len <- cache$len #branch length
  	depth <- cache$depth #branch starting time
  	children <- cache$children #children of each node
  	order <- cache$order[-length(cache$order)] #parent of each node
  	root <- cache$root #root label
	
  	n <- length(len)
  	lq <- rep(0, n)
 	n.tip <- cache$n.tip

  	y <- cache$y #initial tip states
  	branch.init <- branch.base <- vector("list", n)
  	
	#start from tip
    for ( x in 1:n.tip ) {
    	branch.init[[x]] <- y[[x]]
    	ans <- branches(branch.init[[x]], len[x], pars, 0)
    	lq[x] <- ans[[1]]
    	branch.base[[x]] <- c(ans[[2]][1:2],ans[[2]][2]-ans[[2]][1],ans[[2]][3:5])
    }
    
	for ( i in order ) {
		y.in <- initial.conditions(branch.base[children[i,]], pars)
        if ( !all(is.finite(y.in)) )
      		stop("Bad initial conditions: calculation failure along branches?")
      	branch.init[[i]] <- y.in
    	ans <- branches(branch.init[[i]], len[i], pars, depth[i])
    	lq[i] <- ans[[1]]
    	branch.base[[i]] <- ans[[2]]
  	}
	
	y.in <- initial.conditions(branch.base[children[root,]], pars)
   	branch.init[[root]] <- y.in
   	ans <- branches(branch.init[[root]], len[root], pars, depth[root])
    lq[root] <- ans[[1]]
    branch.base[[root]] <- ans[[2]]
    	
	list(init=branch.init, base=branch.base, lq=lq, vals=branch.base[[root]], er.root=branch.base[[children[root,1]]][2])
}   	

make.branches.comp.sp<- function(branches, eps=0) {
    function(y, len, pars, t0, star) {
      ret <- branches(y, len, pars, t0)
      comp.idx <- c(1:length(y))[-c(1,2)]
      if ( all(ret[comp.idx] >= eps) && all(ret[comp.idx] <=1) ) {
      	if ( any(ret[comp.idx] <= 10^-5) && all(ret[comp.idx] > eps)) {
      		q <- sum(ret[comp.idx])
        	ret[comp.idx] <- ret[comp.idx] / q
        	lq <- log(q)
        } else {
        	lq <- 0
        }
        list(lq, ret)
      } else {
        ti <- len[length(len)]/2
        len1 <- c(len[len <= ti], ti)
        len2 <- len[len > ti] - ti
        n1 <- length(len1)

        ret1 <- Recall(y, len1, pars, t0)
        ret2 <- Recall(ret1[[2]][,n1], len2, pars, t0 + ti)
        ret2[[1]] <- ret2[[1]] + ret1[[1]][n1]

        list(c(ret1[[1]][-n1], ret2[[1]]),
             cbind(ret1[[2]][,-n1], ret2[[2]]))
      }
    }
}

make.info.prosse.sp <- function (phy) {
	list(name="prosse_sp",
	     derivs=derivs.prosse_sp,
		 np=3L,
		 ny=6L,
		 k=2L,
		 idx.e=1:2,
		 idx.d=3:6,
		 argnames=default.argnames.prosse(),
		 phy=phy,
		 ml.default="subplex",
		 mcmc.lowerzero=TRUE) 
}

derivs.prosse_sp <- function (t,y,pars) {
			Ei <- y[1]
			Er <- y[2]
			D <- y[3]
			Dr <- y[4]
			Di <- y[5]
			Dir <- y[6]
			b <- pars[1]
			mu <- pars[2]
			lambda <- pars[3]
			c(-(b+mu)*Ei + b*Ei*Ei + mu,
			  -(b+mu+lambda)*Er + b*Er*Er + mu + lambda*Ei,
			  -(b+mu+lambda)*D + b*D*Er,
			  -(b+mu)*Dr + 2*b*Dr*Ei + lambda*D,
			  -(b+mu+lambda)*Di + 2*b*Di*Er + lambda*(D+Dr),
			  -(b+mu+lambda)*Dir + 2*b*Dir*Ei + lambda*(D+Dr))
}
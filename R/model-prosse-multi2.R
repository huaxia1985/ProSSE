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
make.prosse.multi2 <- function (tree, traits, states, states.sd, lambda, branch.init.old=NULL, branch.base.old=NULL, branch.star.old=NULL, lq.old=NULL, Et.old=NULL, change.tip=NULL, control=NULL) {
	cache <- make.cache.prosse.multi(tree, traits, states, states.sd, lambda, control)
    all.branches <- make.all.branches.prosse.multi2(cache, branch.init.old, branch.base.old, branch.star.old, lq.old, Et.old, change.tip, cache$control)
	rootfunc <- make.rootfunc.prosse.multi2(cache)
	f.pars <- make.pars.prosse.multi(cache)
    ll <- function (pars, condition.surv=TRUE, root=ROOT.OBS, root.f=NULL) {
		pars2 <- f.pars(pars)
		ans <- all.branches(pars2)
        rootfunc(ans, pars2, condition.surv, root, root.f)
	}
	class(ll) <- c("prosse","dtlik","function")
	ll
}

make.rootfunc.prosse.multi2 <- function(cache) {
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

    list(loglik=loglik,branch.init=res$init,branch.base=res$base,branch.star=res$star,lq=res$lq, Et=res$Et)
  }
}

######################################################################
## Extra core stuff:
#this is similar to make.all.branches.dtlik
make.all.branches.prosse.multi2 <- function (cache, branch.init.old=NULL, branch.base.old=NULL, branch.star.old=NULL, lq.old=NULL, Et.old=NULL, change.tip=NULL, control) {
	branches <- make.branches.prosse.multi(cache,control)
	initial.conditions <- make.initial.conditions.prosse.multi(control)
	function (pars) {
		cache$y <- initial.tip.prosse.multi(cache, control, pars[[1]]$x)
        all.branches.matrix.prosse.multi2(pars,cache,initial.conditions,branches,branch.init.old,branch.base.old,branch.star.old,lq.old,Et.old,change.tip)
	}
}

all.branches.matrix.prosse.multi2 <- function (pars,cache,initial.conditions,branches, branch.init.old=NULL, branch.base.old=NULL, branch.star.old=NULL, lq.old=NULL, Et.old=NULL, change.tip=NULL) {
	len <- cache$len #branch length
  	depth <- cache$depth #branch starting time
  	children <- cache$children #children of each node
  	order <- cache$order #parent of each node
  	root <- cache$root #root label
	group <- cache$group #list of clades
	stars <- cache$stars #indicator for which ODE to use
	h <- cache$control$h
	
  	n <- length(len)
 	
 	n.tip <- cache$n.tip

  	y <- cache$y #initial tip states
  	
  	if (is.null(change.tip)) {
  	lq <- rep(0, n)  #this is the log of DI+DR after the calculation on each branch. It is used to standardize DI and DR to avoid rounding error due to extremely small likelihood value
  	branch.init <- branch.base <- vector("list", n)
	Et <- rep(0,n)
    branch.star <- rep(0,n)
	
	#start from each clade
    for ( x in 1:length(y) ) {
        idx <- y[[x]]$target #tip index in each clade
        if (length(idx)==1) { #for clade with single tip
        	branch.init[idx] <- list(y[[x]]$y)
        	ans <- branches(y[[x]]$y, y[[x]]$t, pars, 0, star=2)
            branch.star[idx] <- 2
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
                branch.star[idx[j]] <- 1
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
                branch.star[j] <- 1
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
                branch.star[j] <- 2
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
                    branch.star[j] <- 2
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
        branch.star[i] <- 2
        lq[i] <- y.in[1] + ans[[1]]
    	Et[i] <- ans[[2]]
    	branch.base[[i]] <- ans[[3]]
    	}
  	}
  	} else {
  		branch.init <- branch.init.old
  		branch.base <- branch.base.old
        branch.star <- branch.star.old
  		lq <- lq.old
  		Et <- Et.old
	
	#start from each clade
    for ( x in 1:length(y) ) {
        idx <- y[[x]]$target #tip index in each clade
        if (length(idx)==1) {
            branch.star[idx] <- 2
            if (!identical(branch.init[idx],list(y[[x]]$y))||!identical(branch.star[idx],branch.star.old[idx])) { #for clade with single tip
        	branch.init[idx] <- list(y[[x]]$y)
        	ans <- branches(y[[x]]$y, y[[x]]$t, pars, 0, star=2)
        	lq[idx] <- ans[[1]]
        	Et[idx] <- ans[[2]]
        	branch.base[[idx]] <- ans[[3]]
        	}
        } else if (length(idx)>1) { #for clade more than one tip
        	idx2 <- group[[x]]
        	idx2 <- idx2[idx2>n.tip]
        	#calculate along tip branch
                for (j in c(1:length(idx))) {
                    branch.star[idx[j]] <- 1
                    if (!identical(branch.init[idx[j]],y[[x]]$y[j])||!identical(branch.star[idx[j]],branch.star.old[idx[j]])) {
                    branch.init[idx[j]] <- y[[x]]$y[j]
                    ans <- branches(y[[x]]$y[[j]], y[[x]]$t[j], pars, 0, star=1) #use the ODE for DI*, so star=1
                    lq[idx[j]] <- ans[[1]]
                    Et[idx[j]] <- ans[[2]]
                    branch.base[[idx[j]]] <- ans[[3]]
                }
                }
        	#calculate along internal branch
         	for (j in idx2[-length(idx2)]) {
         		tmp <- children[j,!is.element(children[j,],group[[x]])] 
         		if (length(tmp)>0) {
         			branch.base[[tmp]][,(h+1):(2*h)] <- 0
         		} #find children that leads to another species and set its DR=0
                branch.star[j] <- 1
         		if (!identical(branch.base[children[j,]],branch.base.old[children[j,]])||!identical(branch.star[j],branch.star.old[j])) {
        		y.in <- initial.conditions(branch.base[children[j,]], pars, depth[j]) #calculate the boundary condition of the branch
        		if ( !all(is.finite(y.in)) )
      				stop("Bad initial conditions: calculation failure along branches?")
    			branch.init[[j]] <- y.in[-1]
    			ans <- branches(c(Et[children[j,1]], y.in[-1]), len[j], pars, depth[j], star=1)
                lq[j] <- y.in[1] + ans[[1]]
    			Et[j] <- ans[[2]]
    			branch.base[[j]] <- ans[[3]]
    			}
        	}
        	#calculate the boundary conditions for the basal branch
        	j <- idx2[length(idx2)]
        	tmp <- children[j,!is.element(children[j,],group[[x]])] 
         	if (length(tmp)>0) {
         		branch.base[[tmp]][,(h+1):(2*h)] <- 0
         	} #find children that leads to another species and set its DR=0
            branch.star[j] <- 2
         	if (!identical(branch.base[children[j,]],branch.base.old[children[j,]])||!identical(branch.star[j],branch.star.old[j])) {
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
        	}
        } else { #if the group is not a clade, but the branch linking a clade to another clade that it is nested in, then we need to calculate this branch first.
        	for (j in group[[x]]) {
                branch.star[j] <- 2
        		if (!identical(branch.base[children[j,]],branch.base.old[children[j,]])||!identical(branch.star[j],branch.star.old[j])) {
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
    }

    for ( i in order ) {
    	if (i!=root) {
            branch.star[i] <- 2
    	if (!identical(branch.base[children[i,]],branch.base.old[children[i,]])||!identical(branch.star[i],branch.star.old[i])) {
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
  	}
  	}
 	
  	y.in <- initial.conditions(branch.base[children[root,]], pars, depth[root])
  	branch.init[[root]] <- y.in[-1]
  	ans <- branches(c(Et[children[root,1]], y.in[-1]), len[root], pars, depth[root], star=2)
    branch.star[root] <- 2
    lq[root] <- y.in[1] + ans[[1]]
    Et[root] <- ans[[2]]
    branch.base[[root]] <- ans[[3]]
  	
    list(init=branch.init, base=branch.base, star=branch.star, lq=lq, Et=Et, vals=branch.base[[root]], e.root=Et[root])
}

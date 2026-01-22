make.prosse.trait.lambda <- function (tree,traits,types) {
	cache <- make.cache.prosse.trait.lambda(tree,traits,types)
	f.pars <- make.pars.prosse.trait.lambda(cache)
	all.branches <- make.all.branches.prosse.trait.lambda(cache)
	rootfunc <- function(res, pars, t=max(cache$depth),k=cache$k, condition.surv, root, root.f) {
  		h <- prod(k)
  		vals <- res$vals
    	lq <- res$lq
   		d.root <- res$vals[(h+2):(2*h+1)]
			
		if (condition.surv) {
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
	
    	root.p <- diversitree:::root.p.calc(d.root, pars, root)
		loglik <- log(sum(root.p * d.root)) + sum(lq)

    	loglik
  	}
	ll <- function (pars, condition.surv=T,root=ROOT.OBS,root.f=NULL) {
		pars2 <- f.pars(pars)
		ans <- all.branches(pars2)
		rootfunc(ans, pars, max(cache$depth), cache$k, condition.surv, root, root.f)
	}
	class(ll) <- c("prosse","dtlik","function")
	ll
}

## 2: info
make.info.prosse.trait.lambda <- function (k, phy) {
	list(name="prosse_trait_lambda",
		 ## Parameters:
		 np=NA,
		 argnames=default.argnames.prosse.trait.lambda(k),
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

default.argnames.prosse.trait.lambda <- function(k) {
  c("b","mu","lambda.m","lambda.c",
    unlist(lapply(1:length(k), function (j) sprintf("q%s%s.%s", rep(as.character(1:k[j]), each=k[j]-1), unlist(lapply(1:k[j], function(i) as.character(1:k[j])[-i])),j))),
    unlist(lapply(1:length(k), function (j) sprintf("p%s%s.%s", rep(as.character(1:k[j]), each=k[j]-1), unlist(lapply(1:k[j], function(i) as.character(1:k[j])[-i])),j)))
   )
}


## 3: make.cache (& initial conditions)
make.cache.prosse.trait.lambda <- function (tree, traits, types) {
  	#tree <- check.tree(tree)
	cache <- make.cache.tree.prosse.lambda(tree,types)
    #trait parameters
    k <- sapply(1:ncol(traits),function (i) length(unique(traits[,i])))
    cache$traits <- as.matrix(traits)
    cache$k <- k
    cache$h <- prod(k)
    n.q <- sum(k*(k-1))
    n.args   <- 4 + n.q + n.q
    args <- list(b = 1,
                 mu = 2,
                 lambda.m = 3,
                 lambda.c = 4,
                 q = 4 + seq_len(n.q),
                 p = 4 + n.q + seq_len(n.q)
                 )
    cache$args <- args
    #info
   	cache$info <- make.info.prosse.trait.lambda(k, tree)
	cache
}

## 4: make.pars
make.pars.prosse.trait.lambda <- function (cache) {
	args <- cache$args
	k <- cache$k
	h <- cache$h
	function (pars) {
		b <- pars[args$b]
   		mu <- pars[args$mu]
   		lambda.m <- pars[args$lambda.m]
   		lambda.c <- pars[args$lambda.c]
		q <- pars[args$q]
		p <- pars[args$p]
		#calculate q_r matrix
    	nam <- expand.grid(lapply(k,function (i) 1:i))
    	q.r <- diag(dim(nam)[1])
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
    	q.m <- diag(dim(nam)[1])
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
    	names(pars) <- NULL
		list(b=b,mu=mu,lambda.m=lambda.m,lambda.c=lambda.c,q.r=q.r,q.m=q.m)
	}
}

make.all.branches.prosse.trait.lambda <- function (cache) {
	h <- cache$h
	branches <- function(y, len, pars, t0, star, type) {
		
   		#calculate e^(Q*dt)
    	Qr <- as.matrix(Matrix::expm(pars$q.r*len))
		Q <- as.matrix(Matrix::expm(pars$q.m*len))

		b <- pars$b
		mu <- pars$mu
		lambda <- ifelse(type==0,pars$lambda.m,pars$lambda.c)
		r <- b-mu
		z <- exp(len * r)
		z2 <- exp(- len * lambda)
		e0 <- y[1]
  		y[1] <- (mu + z*(e0 - 1)*mu - b*e0) / (mu + z*(e0 - 1)*b - b*e0)
  		A <- (z * r * r)/(z * b - mu + (1-z)*b*e0)^2
  		if (star==1) {
  			y[2:(h+1)] <- Q %*% (A * y[2:(h+1)] * z2)
  			y[(h+2):(2*h+1)] <- Q %*% (A * y[(h+2):(2*h+1)] * z2)
  		} else if (star==2) {
  			if (type==0) {
  				dr0 <- A * y[(h+2):(2*h+1)]
    			y[(h+2):(2*h+1)] <- Qr %*% (A * y[(h+2):(2*h+1)])
    			d0 <- A * y[2:(h+1)] * z2 - dr0 * z2
    			y[2:(h+1)] <- y[(h+2):(2*h+1)] + Q %*% d0
    		} else {
    			y[2:(h+1)] <- Q %*% (A * y[2:(h+1)] * z2) - A * (z2 - 1) * y[(h+2):(2*h+1)]
    			y[(h+2):(2*h+1)] <- Q %*% (A * y[(h+2):(2*h+1)])
    		}
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
		cache$y <- initial.tip.prosse.trait(cache)
		all.branches.matrix.prosse.trait.lambda(pars,cache,initial.conditions=initial.conditions.prosse.trait,branches)
	}
}

all.branches.matrix.prosse.trait.lambda <- function (pars,cache,initial.conditions,branches) {
	len <- cache$len #branch length
  	depth <- cache$depth #branch starting time
  	children <- cache$children #children of each node
  	order <- cache$order #parent of each node
  	root <- cache$root #root label
	group <- cache$group #list of clades
	types <- cache$types #indicator for which lambda to use
	h <- cache$h
	
  	n <- length(len)
  	lq <- rep(0, n)  #this is the log of DI+DR after the calculation on each branch. It is used to standardize DI and DR to avoid rounding error due to extremely small likelihood value
 	n.tip <- cache$n.tip

  	y <- cache$y #initial tip states
  	branch.init <- branch.base <- vector("list", n)
	
	#start from each clade
    for ( x in 1:length(y) ) {
        idx <- y[[x]]$target #tip index in each clade
        if (length(idx)==1) { #for clade with single tip
        	branch.init[idx] <- list(y[[x]]$y)
        	ans <- branches(y[[x]]$y, y[[x]]$t, pars, 0, star=2, type=types[idx])
        	lq[idx] <- ans[1]
        	branch.base[[idx]] <- ans[-1]
        } else if (length(idx)>1){ #for clade more than one tip
        	idx2 <- group[[x]]
        	idx2 <- idx2[idx2>n.tip]
        	#calculate along tip branch
        	for (j in 1:length(idx)) { 
        		branch.init[idx[j]] <- y[[x]]$y[j]
        		ans <- branches(y[[x]]$y[[j]], y[[x]]$t[j], pars, 0, star=1, type=types[idx[j]]) #use the ODE for DI*, so star=1
         		lq[idx[j]] <- ans[1]
        		branch.base[[idx[j]]] <- ans[-1]
        	}
        	#calculate along internal branch
         	for (j in idx2[-length(idx2)]) {
         		tmp <- children[j,!is.element(children[j,],group[[x]])] 
         		if (length(tmp)>0) {
         			if (is.matrix(branch.base[[tmp]])) {
         				branch.base[[tmp]][,(h+2):(2*h+1)] <- 0
         			} else {
         				branch.base[[tmp]][(h+2):(2*h+1)] <- 0
         			}
         		} #find children that leads to another species and set its DR=0
        		y.in <- initial.conditions(branch.base[children[j,]], pars$b, depth[j], h) #calculate the boundary condition of the branch
        		if ( !all(is.finite(y.in)) )
      				stop("Bad initial conditions: calculation failure along branches?")
    			branch.init[[j]] <- y.in
    			ans <- branches(y.in, len[j], pars, depth[j], star=1, type=types[j])
    			lq[j] <- ans[1]
    			branch.base[[j]] <- ans[-1]
        	}
        	#calculate the boundary conditions for the basal branch
        	j <- idx2[length(idx2)]
        	tmp <- children[j,!is.element(children[j,],group[[x]])] 
         	if (length(tmp)>0) {
         		if (is.matrix(branch.base[[tmp]])) {
         				branch.base[[tmp]][,(h+2):(2*h+1)] <- 0
         			} else {
         				branch.base[[tmp]][(h+2):(2*h+1)] <- 0
         			}
         	} #find children that leads to another species and set its DR=0
        	y.in <- initial.conditions(branch.base[children[j,]], pars$b, depth[j], h)
        	if ( !all(is.finite(y.in)) )
      			stop("Bad initial conditions: calculation failure along branches?")
    		branch.init[[j]] <- y.in
    		#calculate along the basal branch
    		if (j!=root) {
        		ans <- branches(y.in, len[j], pars, depth[j], star=2, type=types[j]) #use the ODE for DI, so star=2
        		lq[j] <- ans[1]
        		branch.base[[j]] <- ans[-1]
        	}
        } else { #if the group is not a clade, but the branch linking a clade to another clade that it is nested in, then we need to calculate this branch first.
        	for (j in group[[x]]) {
        		y.in <- initial.conditions(branch.base[children[j,]], pars$b, depth[j], h)
        		if ( !all(is.finite(y.in)) )
      				stop("Bad initial conditions: calculation failure along branches?")
    			branch.init[[j]] <- y.in
    			if (j!=root) {
    				ans <- branches(y.in, len[j], pars, depth[j], star=2, type=types[j])
    				lq[j] <- ans[1]
        			branch.base[[j]] <- ans[-1]
        		}
        	}
        }
    }

    for ( i in order ) {
    	if (i!=root) {
    	y.in <- initial.conditions(branch.base[children[i,]], pars$b, depth[i], h)
    	if ( !all(is.finite(y.in)) )
      		stop("Bad initial conditions: calculation failure along branches?")
    	branch.init[[i]] <- y.in
    	ans <- branches(y.in, len[i], pars, depth[i], star=2, type=types[i])
    	lq[i] <- ans[1]
    	branch.base[[i]] <- ans[-1]
    	}
  	}
 	
  	y.in <- initial.conditions(branch.base[children[root,]], pars$b, depth[root], h)
  	branch.init[[root]] <- y.in
  	ans <- branches(y.in, len[root], pars, depth[root], star=2, type=0)
    lq[root] <- ans[1]
    branch.base[[root]] <- ans[-1]
  	
  list(init=branch.init, base=branch.base, lq=lq, vals=branch.base[[root]])
}

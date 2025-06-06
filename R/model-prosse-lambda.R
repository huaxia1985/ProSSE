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
make.prosse.lambda <- function (tree,types) {
	#types: 0=morpho species, 1=species complex 1, 2=species.complex 2, ....
	cache <- make.cache.prosse.lambda(tree,types)
	all.branches <- make.all.branches.prosse.lambda(cache)
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
		ans <- all.branches(pars)
		rootfunc(ans, pars, max(cache$depth), condition.surv)
	}
	class(ll) <- c("prosse","dtlik","function")
	ll
}

## 2: info
make.info.prosse.lambda <- function (phy) {
	list(name="prosse_lambda",
		 ## Parameters:
		 np=4L,
		 argnames=default.argnames.prosse.lambda(),
		 ## Variables:
		 ny=4L,
		 k=NA,
		 idx.e=NA,
		 idx.d=NA,
		 ## Phylogeny:
		 phy=phy,
		 ## Inference:
		 ml.default="subplex",
		 mcmc.lowerzero=TRUE)
	}

default.argnames.prosse.lambda <- function() {
	c("b","mu","lambda.m","lambda.c")
}

## 3: make.cache (& initial conditions)
make.cache.prosse.lambda <- function (tree,types) {
  	#tree <- check.tree(tree)
	cache <- make.cache.tree.prosse.lambda(tree,types)
	cache$info <- make.info.prosse.lambda(tree)
	cache
}

######################################################################
## Extra core stuff:
#this is similar to make.all.branches.dtlik
make.all.branches.prosse.lambda <- function (cache) {
	branches <- function(y, len, pars, t0, star, type) {
		b <- pars[1]
		mu <- pars[2]
		lambda <- ifelse(type==0,pars[3],pars[4])
		r <- b-mu
		z <- exp(len * r)
		z2 <- exp(- len * lambda)
		e0 <- y[1]
  		dr0 <- y[3]
  		y[1] <- (mu + z*(e0 - 1)*mu - b*e0) / (mu + z*(e0 - 1)*b - b*e0)
  		A <- (z * r * r)/(z * b - mu + (1-z)*b*e0)^2
  		if (star==1) {
  			y[-1] <- A * y[-1] * z2
  		} else if (star==2) {
  			y[3] <- A * y[3]
  			y[2] <- y[3] + (y[2]-dr0) * A * z2
  		}
  		if (any(y[-1]<=10^-5 && y[-1]>0)) {
  			lq <- sum(y[-1])
  			y[-1] <- y[-1] / lq
  			lq <- log(lq)
  		} else {
  			lq <- 0
  		}
  		c(lq, y) 		
		}
	function (pars) {
		cache$y <- initial.tip.prosse(cache)
		all.branches.matrix.prosse.lambda(pars,cache,initial.conditions=initial.conditions.prosse,branches)
	}
}

#this is similar to make.branches.matrix
#the core calcuation function 
all.branches.matrix.prosse.lambda <- function (pars,cache,initial.conditions,branches) {
	len <- cache$len #branch length
  	depth <- cache$depth #branch starting time
  	children <- cache$children #children of each node
  	order <- cache$order #parent of each node
  	root <- cache$root #root label
	group <- cache$group #list of clades
	types <- cache$types #indicator for which lambda to use
	
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
         				branch.base[[tmp]][,3] <- 0
         			} else {
         				branch.base[[tmp]][3] <- 0
         			}
         		} #find children that leads to another species and set its DR=0
        		y.in <- initial.conditions(branch.base[children[j,]], pars, depth[j]) #calculate the boundary condition of the branch
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
         				branch.base[[tmp]][,3] <- 0
         			} else {
         				branch.base[[tmp]][3] <- 0
         			}
         	} #find children that leads to another species and set its DR=0
        	y.in <- initial.conditions(branch.base[children[j,]], pars, depth[j])
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
        		y.in <- initial.conditions(branch.base[children[j,]], pars, depth[j])
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
    	y.in <- initial.conditions(branch.base[children[i,]], pars, depth[i])
    	if ( !all(is.finite(y.in)) )
      		stop("Bad initial conditions: calculation failure along branches?")
    	branch.init[[i]] <- y.in
    	ans <- branches(y.in, len[i], pars, depth[i], star=2, type=types[i])
    	lq[i] <- ans[1]
    	branch.base[[i]] <- ans[-1]
    	}
  	}
 	
  	y.in <- initial.conditions(branch.base[children[root,]], pars, depth[root])
  	branch.init[[root]] <- y.in
  	ans <- branches(y.in, len[root], pars, depth[root], star=2, type=0)
    lq[root] <- ans[1]
    branch.base[[root]] <- ans[-1]
  	
  list(init=branch.init, base=branch.base, lq=lq, vals=branch.base[[root]])
}

make.cache.tree.prosse.lambda <- function(tree,types) {
  if (inherits(tree, "phylo"))
    class(tree) <- "phylo"
  edge <- tree$edge
  edge.length <- tree$edge.length
  idx <- seq_len(max(edge))
  n.tip <- length(tree$tip.label)
  tips <- seq_len(n.tip)
  root <- n.tip + 1
  
  is.tip <- idx <= n.tip

  children <- diversitree:::get.children(edge, n.tip)
  parent <- edge[match(idx, edge[,2]),1]
  order <- order.tmp <- diversitree:::get.ordering(children,is.tip,root)
  
  species <- tree$species
  
  species.label <- unique(species)
  n.species <- length(species.label)
  tip.species <- sapply(species.label,function(i) is.element(species,i))
  n.tip.species <- colSums(tip.species)
  
  #for species as 1 tip
  group <- lapply(which(n.tip.species==1),function (i) which(tip.species[,i]))
 
  #for species as monophyletic clade with more than 1 tips
  clade <- which(n.tip.species>1)
  if (length(clade)>0) {
  mono <- sapply(clade,function (i) is.monophyletic(tree,which(tip.species[,i])))
  if (length(mono)>0) {
  group.tmp <- lapply(clade[mono],function (i) edge[max(which(is.element(edge[,2],which(tip.species[,i]))))-c(0:(n.tip.species[i]*2-2)),2])
  group <- c(group,group.tmp)
  }
  
  #for species as paraphyletic clade
  anc <- vector("list", max(order))
  for (i in c(rev(order[-length(order)]),tips)) {anc[[i]] <- c(parent[i],anc[[parent[i]]])}
  exc <- lapply(clade[!mono],function (i) Reduce(intersect,anc[which(tip.species[,i])]))
  if (length(exc)>0) {
  nodes <- lapply(1:length(exc), function (i) unlist(anc[which(tip.species[,clade[!mono][i]])]))
  group.tmp <- lapply(1:length(exc), function (i) unique(c(which(tip.species[,clade[!mono][i]]),sort(nodes[[i]][!is.element(nodes[[i]],exc[[i]][-1])],decreasing=T))))
  group.tmp2 <- unlist(group.tmp)
  group.root <- sapply(exc,max,na.rm=T)
  }
  names(group.tmp) <- names(exc)
  
  if (length(exc)>0) {
  #order calculation from 1-tip species, monophyletic species, any branch that connects these species to the root of a paraphyletic species if nested.
  tmp <- lapply(group,function (i) anc[[i[length(i)]]])
  tmp2 <- lapply(tmp,function (i) which(is.element(i,group.tmp2)))
  tmp3 <- which(sapply(tmp2,length)>0)
  if (length(tmp3)>0) {
    tmp <- unlist(lapply(tmp3,function (i) tmp[[i]][0:(tmp2[[i]][1]-1)]))
    tmp <- order[is.element(order,tmp)]
  } else {
  	tmp <- NULL
  }
  
  #order calculation from paraphyletic species, any branch that connects these species to the root of a paraphyletic species if nested.
  tmp_para <- lapply(group.tmp,function (i) anc[[i[length(i)]]])
  tmp2 <- lapply(1:length(tmp_para),function (i) which(is.element(tmp_para[[i]],unlist(group.tmp[-i]))))
  tmp3 <- which(sapply(tmp2,length)>0)
  if (length(tmp3)>0) {
    tmp_para <- unlist(lapply(tmp3,function (i) tmp_para[[i]][0:(tmp2[[i]][1]-1)]))
    tmp_para <- order[is.element(order,tmp_para)]
  } else {
  	tmp_para <- NULL
  }
  tmp <- unique(c(tmp,tmp_para))
  if (length(tmp)>0) {
  	tmp <- as.list(tmp) 
  }

  group <- c(group,tmp,group.tmp)
  
  #check if a paraphyletic group contains another paraphyletic group is correct
  if (length(tmp)+(length(group.tmp))>1) {
  	i <- length(group)-length(group.tmp)-length(tmp)+1
  	while (i < length(group)) {
  		x <- 0
        tmp2 <- group[[i]]
        tmp2 <- tmp2[tmp2>n.tip]
        #check if each internal branch has children calculated
        for (j in 1:length(tmp2)) { 
       		if (prod(is.element(children[tmp2[j],],unlist(group[1:i])))==0) {
       			group <- c(group[-i],group[i])
       			x <- 1
       			break
       		}
        }
        if (x ==0) {i <- i+1}
     }
  }

  order <- order[!is.element(order,c(unlist(group),tmp))]
  }
  }
  
  len <- edge.length[match(idx, edge[,2])]

  height <- diversitree:::branching.heights(tree)
  depth <- max(height) - height
  len[root] <- tree$root.depth-max(depth)

  if ( is.ultrametric(tree) )
    depth[tips] <- 0

  #for species complex
  types.tmp <- c(types,numeric(n.tip-1)) #0 for lambda.m, 1 for lambda.c
  for (i in 1:max(types)) {
  	tip <- which(types==i)
  	node <- as.numeric(names(which(table(unlist(anc[tip])) < length(tip))))
  	types.tmp[node] <- i
  }
  
  ans <- list(tip.label=tree$tip.label,
              node.label=tree$node.label,
              len=len,
              children=children,
              parent=parent,
              order=order,
              root=root,
              n.tip=n.tip,
              n.node=tree$Nnode,
              tips=tips,
              height=height,
              depth=depth,
              edge=edge,
              edge.length=edge.length,
              group=group,
              types=types.tmp,
              species=species)
  ans
}


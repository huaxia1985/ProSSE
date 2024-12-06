## MCMC for prosse, using slice sampler for pars and gibbs sampler for node state

mcmc.prosse <- function(lik, ..., tree, species.name, unknown.tip, unknown.list=NULL, traits=NULL, types=NULL, states=NULL, states.sd=NULL, lambda=NULL, control=NULL, x.init, nstepsw=0, nsteps, w, prior=NULL,
                         lower=-Inf, upper=Inf, sampler=sampler.slice, fail.value=-Inf,
                         print.every=1,
                         save.file, save.every=0, save.every.dt=NULL,
                         previous=NULL,previous.tol=1e-4,keep.func=TRUE) {
  if ( is.null(sampler) )
    sampler <- diversitree:::sampler.slice
  
  if ( save.every > 0 || !is.null(save.every.dt) ) {
    if ( missing(save.file) )
      stop("save.file must be given if save.every > 0")
    save.type <- tools::file_ext(save.file)
    save.type <- match.arg(tolower(save.type), c("csv", "rds"))
    save.fun <- switch(save.type,
                       rds=saveRDS,
                       csv=function(x, f) write.csv(x, f, row.names=FALSE))
    save.file.bak <- paste(save.file, ".bak", sep="")
  }

  n.previous <- if ( is.null(previous) ) 0 else nrow(previous)
  
  par.args <- c(argnames(lik), unknown.tip)
  n.par <- length(argnames(lik))
  npar <- n.par + length(unknown.tip)
  
  lower <- diversitree:::check.par.length(lower, n.par)
  upper <- diversitree:::check.par.length(upper, n.par)
  w     <- diversitree:::check.par.length(w,     n.par)
  
  if ( is.null(prior) ) {
       posterior <- protect(function(x) lik(x),
                          fail.value.default=fail.value)
   } else {
       if (is.list(prior)) {
           posterior <- protect(function(x) lik(x)+ sum(sapply(1:length(prior),function (i) prior[[i]](x[i]))),
           fail.value.default=fail.value)
       } else {
           posterior <- protect(function(x) lik(x)+ prior(x),
                          fail.value.default=fail.value)
       }
   }
   
  if ( !is.null(previous) ) {
      if ( !inherits(previous$out, "mcmcsamples") )
        stop("Currently only mcmcsamples objects can be continued")
      if ( n.previous >= nsteps ) {
        warning("Chain already complete")
        return(previous)
      }
    if ( !is.null(x.init) )
      stop("x.init must be NULL if continuing")
    hist.pars <- matrix(NA, ncol=npar, nrow=nsteps)
    hist.prob <- rep(NA, nsteps)
    
    hist.pars[seq_len(n.previous),] <- stats::coef(previous)
    hist.prob[seq_len(n.previous)]  <- previous$p

    x.init <- hist.pars[n.previous,]
    y.prev <- hist.prob[n.previous]
  } else {
  	hist.pars <- matrix(NA, ncol=npar, nrow=nsteps+nstepsw)
    hist.prob <- rep(NA, nsteps+nstepsw)
      if ( is.null(names(x.init)) ) {
        try(colnames(hist.pars) <- names(x.init) <- par.args,
            silent=TRUE)
      } else {
        colnames(hist.pars) <- names(x.init)
      }
    diversitree:::check.bounds(lower, upper, x.init[1:n.par])
    y.init <- posterior(x=x.init[1:n.par], fail.value=fail.value)
  }
  if ( !is.finite(y.init) ||
      (!is.null(fail.value) && y.init == fail.value) )
    stop("Starting point must have finite probability")
  if ( n.previous > 0 && abs(y.prev - y.init) > previous.tol ) {
    msg <- paste("Cannot continue the chain:\n",
                 sprintf("\texpected posterior = %2.7, got %2.7f",
                         previous$p[n.previous], y.init))
    stop(msg)
  }
  
  class.str <- c(sprintf("mcmcsamples.%s", diversitree:::get.info(lik)$name),
                 "mcmcsamples", "data.frame")

  clean.hist <- function(pars, p) {
    out <- data.frame(i=seq_along(p), pars, p)
    class(out) <- class.str
    out
  }

  we.should.print <- diversitree:::make.every.so.often(print.every)
  we.should.save <- diversitree:::make.every.so.often(save.every, save.every.dt)
	
  constraint <- list(...)
  		
  mcmc.loop <- function() {
  	if (nstepsw>0) {
  	for ( i in seq(1, nstepsw, by=1) ) {
      lik.tmp <- y.init
      for (j in unknown.tip) {
      	if (is.null(unknown.list)) {
          species <- unique(c(tree$species,species.name[j]))
         } else {
          species <- unique(c(tree$species[unknown.list[[j]]],species.name[j]))
         }
          nspecies <- length(species)
          tmp <- numeric(nspecies)
          zz <- which(species==tree$species[j])
          tmp[zz] <- lik.tmp
          for (z in c(1:nspecies)[-zz]) {
              tree$species[j] <- species[z]
              if (check.paraphyletic(tree)) {
                  tmp[z] <- fail.value
              } else {
                  lik <- make.prosse(tree, traits, types, states, states.sd, lambda, control)
                  if (!is.null(constraint)) {
                  	lik <- constrain(lik,...)
                  }
                  if ( is.null(prior) ) {
                       posterior <- protect(function(x) lik(x),
                                          fail.value.default=fail.value)
                   } else {
                       posterior <- protect(function(x) lik(x)+ prior(x),
                                          fail.value.default=fail.value)
                   }
                  tmp[z] <- posterior(x=x.init[1:n.par], fail.value=fail.value)
              }
          }
          prob <- exp(tmp-max(tmp))
          prob <- prob/sum(prob)
          x.init[j] <- hist.pars[i,j] <- tree$species[j] <- z <- sample(x=species,size=1,prob=prob)
          lik.tmp <- tmp[which(species==z)]
      }
      y.init <- lik.tmp
      lik <- make.prosse(tree, traits, types, states, states.sd, lambda, control)
      if (!is.null(constraint)) {
             lik <- constrain(lik,...)
      }
      if ( is.null(prior) ) {
           posterior <- protect(function(x) lik(x),
                              fail.value.default=fail.value)
       } else {
           posterior <- protect(function(x) lik(x)+ prior(x),
                              fail.value.default=fail.value)
       }
      tmp <- sampler(posterior, x.init[1:n.par], y.init, w,
                     lower, upper, control)
      x.init[1:n.par] <- hist.pars[i,1:n.par] <- tmp[[1]]
      y.init <- hist.prob[i] <- tmp[[2]]
	 
      if ( we.should.print() )
          cat(sprintf("%d: {%s} -> %2.5f\n", i,
          paste(sprintf("%2.4f", x.init), collapse=", "),
          y.init))
      if ( we.should.save() ) {
        j <- seq_len(i)
        ## Back up the old version to avoid IO errors if the system
        ## fails while saving.
        if ( file.exists(save.file) )
          ok <- file.rename(save.file, save.file.bak)
        ok <- try(save.fun(clean.hist(hist.pars[j,,drop=FALSE], hist.prob[j]),
                           save.file))
        if ( inherits(ok, "try-error") )
          warning("Error while writing progress file (continuing)",
                  immediate.=TRUE)
      }
    }
    tmp <- clean.hist(hist.pars[1:nstepsw,], hist.prob[1:nstepsw])
    w <- diff(sapply(tmp[2:(n.par+1)],quantile,c(0.025,0.975)))
	}
    for ( i in seq(n.previous+nstepsw+1, nsteps+nstepsw, by=1) ) {
        lik.tmp <- y.init
      for (j in unknown.tip) {
      	if (is.null(unknown.list)) {
          species <- unique(c(tree$species,species.name[j]))
         } else {
          species <- unique(c(tree$species[unknown.list[[j]]],species.name[j]))
         }
          nspecies <- length(species)
          tmp <- numeric(nspecies)
          zz <- which(species==tree$species[j])
          tmp[zz] <- lik.tmp
          for (z in c(1:nspecies)[-zz]) {
              tree$species[j] <- species[z]
              if (check.paraphyletic(tree)) {
                  tmp[z] <- fail.value
              } else {	
                  lik <- make.prosse(tree, traits, types, states, states.sd, lambda, control)
                  if (!is.null(constraint)) {
                  	lik <- constrain(lik,...)
                  }
                  if ( is.null(prior) ) {
                      posterior <- protect(function(x) lik(x),
                                         fail.value.default=fail.value)
                  } else {
                      posterior <- protect(function(x) lik(x)+ prior(x),
                                         fail.value.default=fail.value)
                  }
                  tmp[z] <- posterior(x=x.init[1:n.par], fail.value=fail.value)
              }
          }
          prob <- exp(tmp-max(tmp))
          prob <- prob/sum(prob)
          x.init[j] <- hist.pars[i,j] <- tree$species[j] <- z <- sample(x=species,size=1,prob=prob)
          lik.tmp <- tmp[which(species==z)]
      }
      y.init <- lik.tmp
      lik <- make.prosse(tree, traits, types, states, states.sd, lambda, control)
      if (!is.null(constraint)) {
          lik <- constrain(lik,...)
      }
      if ( is.null(prior) ) {
          posterior <- protect(function(x) lik(x),
                             fail.value.default=fail.value)
      } else {
          posterior <- protect(function(x) lik(x)+ prior(x),
                             fail.value.default=fail.value)
      }
      tmp <- sampler(posterior, x.init[1:n.par], y.init, w,
                     lower, upper, control)
      x.init[1:n.par] <- hist.pars[i,1:n.par] <- tmp[[1]]
      y.init <- hist.prob[i] <- tmp[[2]]
	 
      if ( we.should.print() )
          cat(sprintf("%d: {%s} -> %2.5f\n", i,
          paste(sprintf("%2.4f", x.init), collapse=", "),
          y.init))
      if ( we.should.save() ) {
        j <- seq_len(i)
        ## Back up the old version to avoid IO errors if the system
        ## fails while saving.
        if ( file.exists(save.file) )
          ok <- file.rename(save.file, save.file.bak)
        ok <- try(save.fun(clean.hist(hist.pars[j,,drop=FALSE], hist.prob[j]),
                           save.file))
        if ( inherits(ok, "try-error") )
          warning("Error while writing progress file (continuing)",
                  immediate.=TRUE)
      }
    }
    clean.hist(hist.pars, hist.prob)
  }

  mcmc.recover <- function(...) {
    j <- !is.na(hist.prob)
    if ( !any(j) )
      stop("MCMC was stopped before any samples taken")
    hist <- clean.hist(hist.pars[j,], hist.prob[j])
    warning("MCMC was stopped prematurely: ", nrow(hist), "/", nsteps,
            " steps completed.  Truncated chain is being returned.",
            immediate.=TRUE)
    hist
  }

  samples <- tryCatch(mcmc.loop(), interrupt=mcmc.recover)

  if ( save.every > 0 || !is.null(save.every.dt) )
    if ( nrow(samples) == nsteps && file.exists(save.file.bak) )
      file.remove(save.file.bak)
      
    if (keep.func) {
    attr(samples, "func")  <- set.defaults(lik, defaults=list(...))
    attr(samples, "prior") <- prior
  }
  
  samples
}

getMRCA <- function (phy, tip) {
    if (length(tip) < 2) {return(NULL)}
    Ntip <- length(phy$tip.label)
    rootnd <- Ntip + 1L
    pars <- integer(phy$Nnode)
    tnd <- tip
    done_v <- logical(Ntip + phy$Nnode)
    pvec <- integer(Ntip + phy$Nnode)
    pvec[phy$edge[, 2]] <- phy$edge[, 1]
    nd <- tnd[1]
    for (k in 1:phy$Nnode) {
        nd <- pvec[nd]
        pars[k] <- nd
        if (nd == rootnd)
            break
    }
    pars <- pars[1:k]
    mrcind <- integer(max(pars))
    mrcind[pars] <- 1:k
    mrcand <- pars[1]
    for (i in 2:length(tnd)) {
        cnd <- tnd[i]
        done <- done_v[cnd]
        while (!done) {
            done_v[cnd] <- TRUE
            cpar <- pvec[cnd]
            done <- done_v[cpar]
            if (cpar %in% pars) {
                if (cpar == rootnd)
                  return(rootnd)
                if (mrcind[cpar] > mrcind[mrcand])
                  mrcand <- cpar
                done_v[cpar] <- TRUE
                done <- TRUE
            }
            cnd <- cpar
        }
    }
    mrcand
}
keep.tip.fixed <- function (tree, tip) {
    des <- unique(tree$edge[sapply(tip,function (i) which(tree$edge[,2]==i)),1])
    node <- des
    des <- des[des!=(length(tree$tip.label)+1)]
    while(length(des)>0) {
        des <- unique(tree$edge[sapply(des,function (i) which(tree$edge[,2]==i)),1])
        node <- c(node,des)
        des <- des[des!=(length(tree$tip.label)+1)]
    }
    mrca <- getMRCA(tree,tip)
    node[node>=mrca]
}
check.paraphyletic <- function (tree) {
    species <- unique(tree$species)
    nspecies <- length(species)
    clades <- sapply(species,function (i) keep.tip.fixed(tree, which(tree$species==i)))
    tmp <- sapply(1:nspecies, function (i) sapply(clades, function (j) sum(is.element(clades[[i]],j))>0))
    tmp <- sapply(1:nspecies,function(i) sapply(c(1:nspecies)[-i], function (j) tmp[i,j]&&tmp[j,i]))
    sum(tmp)>0
}

make.prior.beta <- function(beta) {
  function(pars)
    sum(dbeta(pars, beta, beta, log=TRUE))
}

## MCMC for prosse.multi, using slice sampler for pars and gibbs sampler for node state

mcmc.prosse.multi <- function(lik, tree, species.name, unknown.tip, unknown.list=NULL, traits, states, states.sd, lambda, constraint, control=NULL, x.init, nstepsw=0, nsteps, w, prior=NULL,
                         sampler=sampler.slice2, fail.value=-Inf,
                         lower=-Inf, upper=Inf, print.every=1,
                         save.file, save.every=0, save.every.dt=NULL,
                         previous=NULL,
                         ...) {                 	
  if ( is.null(sampler) )
    sampler <- sampler.slice2
  
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

  n.previous <- if ( is.null(previous) ) 0 else nrow(previous$out)
  
  par.args <- c(argnames(lik), unknown.tip)
  n.par <- length(argnames(lik))
  npar <- n.par + length(unknown.tip)
  
  lower <- check.par.length(lower, n.par)
  upper <- check.par.length(upper, n.par)
  w     <- check.par.length(w,     n.par)
  
  if ( is.null(prior) ) {
       posterior <- protect2(function(x) lik(x),fail.value.default=fail.value)
   } else {
   	if (is.list(prior)) {
   		posterior <- protect2(function(x) {
                             out <- lik(x)
                             out$loglik <- out$loglik+ sum(sapply(1:length(prior),function (i) prior[[i]](x[i])))
                             out
                          },fail.value.default=fail.value)
   	} else {
   		posterior <- protect2(function(x) {
                             out <- lik(x)
                             out$loglik <- out$loglik+ prior(x)
                             out
                          },fail.value.default=fail.value)
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
    
    hist.pars[seq_len(n.previous),] <- stats::coef(previous$out)
    hist.prob[seq_len(n.previous)]  <- previous$out$p

    x.init <- previous$x.init
    y.init <- previous$y.init
    
    colnames(hist.pars) <- names(x.init)
  } else {
  	hist.pars <- matrix(NA, ncol=npar, nrow=nsteps+nstepsw)
    hist.prob <- rep(NA, nsteps+nstepsw)
      if ( is.null(names(x.init)) ) {
        try(colnames(hist.pars) <- names(x.init) <- par.args,
            silent=TRUE)
      } else {
        colnames(hist.pars) <- names(x.init)
      }
    check.bounds(lower, upper, x.init[1:n.par])
    y.init <- posterior(x=x.init[1:n.par], fail.value=fail.value)
  }

  if ( !is.finite(y.init$loglik) ||
      (!is.null(fail.value) && y.init$loglik == fail.value) )
    stop("Starting point must have finite probability")
  
  class.str <- c(sprintf("mcmcsamples.%s", get.info(lik)$name),
                 "mcmcsamples", "data.frame")

  clean.hist <- function(pars, p) {
    out <- data.frame(i=seq_along(p), pars, p)
    class(out) <- class.str
    out
  }

  we.should.print <- make.every.so.often(print.every)
  we.should.save <- make.every.so.often(save.every, save.every.dt)

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
          tmp <- vector("list",nspecies)
          zz <- which(species==tree$species[j])
          tmp[[zz]] <- lik.tmp
          for (z in c(1:nspecies)[-zz]) {
              tree$species[j] <- species[z]
              if (check.paraphyletic(tree)) {
                  tmp[[z]] <- fail.value
              } else {	
                  lik <- make.prosse.multi2(tree, traits, states, states.sd, lambda, branch.init.old=lik.tmp$branch.init, branch.base.old=lik.tmp$branch.base, branch.star.old=lik.tmp$branch.star, lq.old=lik.tmp$lq, Et.old=lik.tmp$Et, change.tip=which(tree$tip.label==j), control=NULL)
                  lik <- constrain2(lik,constraint)
                  if ( is.null(prior) ) {
       posterior <- protect2(function(x) lik(x),fail.value.default=fail.value)
   } else {
   	if (is.list(prior)) {
   		posterior <- protect2(function(x) {
                             out <- lik(x)
                             out$loglik <- out$loglik+ sum(sapply(1:length(prior),function (i) prior[[i]](x[i])))
                             out
                          },fail.value.default=fail.value)
   	} else {
   		posterior <- protect2(function(x) {
                             out <- lik(x)
                             out$loglik <- out$loglik+ prior(x)
                             out
                          },fail.value.default=fail.value)
    }
   }
                  tmp[[z]] <- posterior(x=x.init[1:n.par], fail.value=fail.value)
              }
          }
          prob <- sapply(1:nspecies,function (i) ifelse(length(tmp[[i]])==1,tmp[[i]],tmp[[i]]$loglik))
          prob <- as.numeric(prob)
          prob <- exp(prob-max(prob))
          prob <- prob/sum(prob)
          print(prob)
          x.init[j] <- hist.pars[i,j] <- tree$species[j] <- z <- sample(x=species,size=1,prob=prob)
          lik.tmp <- tmp[[which(species==z)]]
      }
      y.init <- lik.tmp
      lik <- make.prosse.multi2(tree, traits, states, states.sd, lambda)
      lik <- constrain2(lik,constraint)
      if ( is.null(prior) ) {
       posterior <- protect2(function(x) lik(x),fail.value.default=fail.value)
   } else {
   	if (is.list(prior)) {
   		posterior <- protect2(function(x) {
                             out <- lik(x)
                             out$loglik <- out$loglik+ sum(sapply(1:length(prior),function (i) prior[[i]](x[i])))
                             out
                          },fail.value.default=fail.value)
   	} else {
   		posterior <- protect2(function(x) {
                             out <- lik(x)
                             out$loglik <- out$loglik+ prior(x)
                             out
                          },fail.value.default=fail.value)
    }
   }
      tmp <- sampler(posterior, x.init[1:n.par], y.init, w,
                     lower, upper, control)
      x.init[1:n.par] <- hist.pars[i,1:n.par] <- tmp[[1]]
      y.init <- tmp[[2]]
	  hist.prob[i]  <- y.init$loglik
	 
      if ( we.should.print() )
          cat(sprintf("%d: {%s} -> %2.5f\n", i,
          paste(sprintf("%2.4f", x.init), collapse=", "),
          y.init$loglik))
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
    w <- diff(sapply(tmp[2:9],quantile,c(0.025,0.975)))
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
          tmp <- vector("list",nspecies)
          zz <- which(species==tree$species[j])
          tmp[[zz]] <- lik.tmp
          for (z in c(1:nspecies)[-zz]) {
              tree$species[j] <- species[z]
              if (check.paraphyletic(tree)) {
                  tmp[[z]] <- fail.value
              } else {	
                  lik <- make.prosse.multi2(tree, traits, states, states.sd, lambda, branch.init.old=lik.tmp$branch.init, branch.base.old=lik.tmp$branch.base, branch.star.old=lik.tmp$branch.star, lq.old=lik.tmp$lq, Et.old=lik.tmp$Et, change.tip=which(tree$tip.label==j), control=NULL)
                  lik <- constrain2(lik,constraint)
                  if ( is.null(prior) ) {
       posterior <- protect2(function(x) lik(x),fail.value.default=fail.value)
   } else {
   	if (is.list(prior)) {
   		posterior <- protect2(function(x) {
                             out <- lik(x)
                             out$loglik <- out$loglik+ sum(sapply(1:length(prior),function (i) prior[[i]](x[i])))
                             out
                          },fail.value.default=fail.value)
   	} else {
   		posterior <- protect2(function(x) {
                             out <- lik(x)
                             out$loglik <- out$loglik+ prior(x)
                             out
                          },fail.value.default=fail.value)
    }
   }
                  tmp[[z]] <- posterior(x=x.init[1:n.par], fail.value=fail.value)
              }
          }
          prob <- sapply(1:nspecies,function (i) ifelse(length(tmp[[i]])==1,tmp[[i]],tmp[[i]]$loglik))
          prob <- as.numeric(prob)
          prob <- exp(prob-max(prob))
          prob <- prob/sum(prob)
          x.init[j] <- hist.pars[i,j] <- tree$species[j] <- z <- sample(x=species,size=1,prob=prob)
          lik.tmp <- tmp[[which(species==z)]]
      }
      y.init <- lik.tmp
      lik <- make.prosse.multi2(tree, traits, states, states.sd, lambda)
      lik <- constrain2(lik,constraint)
      if ( is.null(prior) ) {
       posterior <- protect2(function(x) lik(x),fail.value.default=fail.value)
   } else {
   	if (is.list(prior)) {
   		posterior <- protect2(function(x) {
                             out <- lik(x)
                             out$loglik <- out$loglik+ sum(sapply(1:length(prior),function (i) prior[[i]](x[i])))
                             out
                          },fail.value.default=fail.value)
   	} else {
   		posterior <- protect2(function(x) {
                             out <- lik(x)
                             out$loglik <- out$loglik+ prior(x)
                             out
                          },fail.value.default=fail.value)
    }
   }
      tmp <- sampler(posterior, x.init[1:n.par], y.init, w,
                     lower, upper, control)
      x.init[1:n.par] <- hist.pars[i,1:n.par] <- tmp[[1]]
      y.init <- tmp[[2]]
	  hist.prob[i]  <- y.init$loglik
	 
      if ( we.should.print() )
          cat(sprintf("%d: {%s} -> %2.5f\n", i,
          paste(sprintf("%2.4f", x.init), collapse=", "),
          y.init$loglik))
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
    out <- clean.hist(hist.pars, hist.prob)
    list(out=out,x.init=x.init,y.init=y.init)
  }

  mcmc.recover <- function(...) {
    j <- !is.na(hist.prob)
    if ( !any(j) )
      stop("MCMC was stopped before any samples taken")
    out <- clean.hist(hist.pars[j,], hist.prob[j])
    warning("MCMC was stopped prematurely: ", nrow(hist), "/", nsteps,
            " steps completed.  Truncated chain is being returned.",
            immediate.=TRUE)
    list(out=out,x.init=x.init,y.init=y.init)
  }

  samples <- tryCatch(mcmc.loop(), interrupt=mcmc.recover)

  if ( save.every > 0 || !is.null(save.every.dt) )
    if ( nrow(samples$out) == nsteps && file.exists(save.file.bak) )
      file.remove(save.file.bak)
      
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

protect2 <- function(f, fail.value.default=NULL) {
  function(..., fail.value=fail.value.default, finite=TRUE) {
    if ( is.null(fail.value) )
      f(...)
    else {
      ret <- try(f(...), silent=TRUE)
      failed <- (inherits(ret, "try-error") ||
                 (finite && !is.finite(ret$loglik)))
      if ( failed )
      	ret <- list(loglik=fail.value)
      else
        ret
    }
  }
}

make.prior.beta <- function(beta) {
  function(pars)
    sum(dbeta(pars, beta, beta, log=TRUE))
}

constrain2 <- function (f, constraint, formulae = NULL, names = argnames(f), extra = NULL) 
{
    if (is.constrained(f)) {
        formulae <- c(attr(f, "formulae"), formulae)
        f <- attr(f, "func")
    }
    formulae <- c(formulae, constraint)
    names.lhs <- names.rhs <- names
    rels <- list()
    for (formula in formulae) {
        res <- constrain.parse(formula, names.lhs, names.rhs, 
            extra)
        if (attr(res, "lhs.is.target")) {
            i <- try(which(sapply(rels, function(x) identical(x, 
                res[[1]]))), silent = TRUE)
            if (inherits(i, "try-error")) 
                stop(sprintf("Error parsing constraint with %s on lhs", 
                  as.character(res[[1]])))
            rels[i] <- res[[2]]
            lhs.txt <- as.character(res[[1]])
            if (any(sapply(rels, function(x) lhs.txt %in% all.vars(x)))) 
                stop(sprintf("lhs (%s) is in an expression and can't be constrained", 
                  lhs.txt))
        }
        names.lhs <- setdiff(names.lhs, unlist(lapply(res, all.vars)))
        names.rhs <- setdiff(names.rhs, as.character(res[[1]]))
        rels <- c(rels, structure(res[2], names = as.character(res[[1]])))
    }
    i <- match(unique(sapply(rels, as.character)), extra)
    final <- c(extra[sort(i[!is.na(i)])], names.rhs)
    npar <- length(final)
    free <- setdiff(names.rhs, names(rels))
    free.i <- match(free, names)
    free.j <- match(free, final)
    target.i <- match(names(rels), names)
    pars.out <- rep(NA, length(names))
    names(pars.out) <- names
    g <- function(pars, ..., pars.only = FALSE) {
        if (length(pars) != npar) 
            stop(sprintf("Incorrect parameter length: expected %d, got %d", 
                npar, length(pars)))
        pars.out[free.i] <- pars[free.j]
        e <- structure(as.list(pars), names = final)
        pars.out[target.i] <- unlist(lapply(rels, eval, e))
        if (pars.only) 
            pars.out
        else f(pars.out, ...)
    }
    class(g) <- c("constrained", class(f))
    attr(g, "argnames") <- final
    attr(g, "formulae") <- formulae
    attr(g, "extra") <- extra
    attr(g, "func") <- f
    g
}
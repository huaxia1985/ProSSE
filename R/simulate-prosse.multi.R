#pars <- c(b,mu,lambda.x,char.evol,list(rate matrix for each trait transition),list(speciation prob matrix for each trait transition))
#x0 initial value of the state and of each trait
## char.evol <- make.brownian.with.drift(0, 0.01)
tree.prosse.multi <- function(pars, max.taxa=Inf, max.t=Inf, x0, single.lineage=TRUE, verbose=F) {
    info <- make.tree.prosse.multi(pars, max.taxa, max.t, x0, single.lineage, verbose)
  if ( single.lineage )
    info <- info[-1,]
  phy <- me.to.ape.prosse.multi(info, attr(info,"trait.name"))
  phy$root.depth <- attr(info, "t")
  if ( is.null(phy) )
    phy
  else
  prune.prosse.multi(phy)
}

make.tree.prosse.multi <- function(pars, max.taxa=Inf, max.t=Inf, x0, single.lineage=TRUE, verbose=F, k=500) {
  b <- pars[[1]]
  mu <- pars[[2]]
  r <- b + mu
  pr.speciation <- b/(b + mu)
  lambda <- pars[[3]]
  char   <- pars[[4]]
  ntrait <- length(pars[[5]])
  trait.q <- lapply(pars[[5]],function(i) if(!is.matrix(i)) as.matrix(i,nrow=length(i)) else i )
  trait.r <- lapply(trait.q, function(i) rowSums(i))
  trait.p <- lapply(pars[[6]],function(i) if(!is.matrix(i)) as.matrix(i,nrow=length(i)) else i )
  trait.k <- sapply(trait.q, function (i) dim(i)[1])
  trait.name <- sprintf("trait.%s", 1:ntrait)
  
  if ( single.lineage ) {
      if (ntrait==1) {
          info <- data.frame(idx=1, len=1e-8, parent=0, species=1, state=x0[1], trait.1=x0[-1], extinct=FALSE, split=FALSE)
      } else {
          info <- data.frame(idx=1, len=1e-8, parent=0, species=1, state=x0[1], trait=list(t(x0[-1])), extinct=FALSE, split=FALSE)
      }
  } else {
      if (ntrait==1) {
          info <- data.frame(idx=1:2, len=1e-8, parent=0, species=1, state=x0[1], trait.1=x0[-1], extinct=FALSE, split=FALSE)
      } else {
          info <- data.frame(idx=1:2, len=1e-8, parent=0, species=1, state=x0[1], trait=list(t(x0[-1])), extinct=FALSE, split=FALSE)
      }
  }
  
  lineages <- which(!info$extinct & !info$split)
  n.taxa <- length(lineages)
  t <- 0
  t.left <- max.t
  n.species <- 1
  
  while ( n.taxa <= max.taxa && n.taxa > 0 && t.left > 0 ) {
  	r.n <- r * n.taxa
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt

    if ( t > max.t ) {
      dt <- dt - (t - max.t)
      info$len[lineages] <- info$len[lineages] + dt
      t <- max.t
      
      ##update state and species identity on lineages during dt
      x <- run.until.change.prosse.multi(lineages, info, k, lambda, char, trait.q, trait.p, ntrait, trait.k, trait.r, trait.name, dt, n.species)
      
      break
    }

    info$len[lineages] <- info$len[lineages] + dt
    
    ##update state and species identity on lineages during dt
    x <- run.until.change.prosse.multi(lineages, info, k, lambda, char, trait.q, trait.p, ntrait, trait.k, trait.r, trait.name, dt, n.species)
    info <- x[[1]]
    n.species <- x[[2]]

    ## Pick a lineage for speciation/extinction event to happen to:
    i <- sample(n.taxa, 1)
    lineage <- lineages[i]

    if ( runif(1) < pr.speciation ) {
      ## Speciating:
      if ( n.taxa == max.taxa ) {
        break
      } else {
        info <- speciate.prosse.multi(info, lineage, trait.name)
        lineages <- c(lineages[-i], c(-1,0) + nrow(info))
        n.taxa <- n.taxa + 1
      }
    } else {
        info$extinct[lineage] <- TRUE
        lineages <- lineages[-i]
        n.taxa <- n.taxa - 1
    }
  }

  attr(info, "t") <- t
  attr(info, "trait.name") <- trait.name
  info
}
  
run.until.change.prosse.multi <- function(lineages, info, k, lambda, char, trait.q, trait.p, ntrait, trait.k, trait.r, trait.name, max.t, n.species) {
    t <- 0
    n.extant <- length(lineages)
    p.change <- 1/k
    repeat {
      state <- info$state[lineages]
      lx <- lambda(state)
      n.state <- lapply(1:ntrait,function (i) sapply(1:trait.k[[i]], function (j) sum(info[lineages,trait.name[i]]==j)))
      r.n <- lapply(1:ntrait,function (i) trait.r[[i]] * n.state[[i]])
      r.tot <- sum(unlist(r.n))
      
      prob <- c(sum(lx), r.tot)
      r <- sum(prob)
      type <- sample(2,1,FALSE,prob)
      dt <- 1/(r*k)
      t <- t + dt
            
      if ( t > max.t ) {
      dt <- dt - (t - max.t)
      }
      
      if ( runif(1) < p.change ) {
        if ( type == 1 ) { # speciation completion due to incompatibility
            i <- sample(n.extant, 1, prob=lx)
            info$species[lineages[i]] <- n.species <- n.species + 1
        } else {
            z <- sample(ntrait,1,FALSE,sapply(r.n,function (i) sum(i)))
                    state <- sample(trait.k[[z]], 1, FALSE, r.n[[z]])
                    lineages.state <- which(!info$extinct & !info$split & info[,trait.name[z]]==state)
                    i <- sample(n.state[[z]][state], 1)
                    if (trait.k[[z]]>2) {
                    	state.new <- sample(c(1:trait.k[[z]])[-state],1,FALSE,trait.q[[z]][state,])
                    } else {
                    	state.new <- c(1,2)[-state]
                    }
                    info[lineages.state[i],trait.name[z]] <- state.new
                    if ( runif(1) < trait.p[[z]][state,state.new] ) {#transition co-occur with speciation completion
                        info[lineages.state[i],"species"] <- n.species <- n.species + 1
                    }
        }
      }

      info$state[lineages] <- char(info$state[lineages], dt)
		
      if ( t >= max.t )
        break
    }

    list(info, n.species)
}
                               
speciate.prosse.multi <- function(info, i, trait.name) {
        j <- 1:2 + nrow(info)
        info[j,"idx"] <- j
        info[j[1],colnames(info)[-1]] <- info[j[2],colnames(info)[-1]] <- c(0,i,info$species[i],info$state[i],info[i,trait.name],FALSE,FALSE)
        info$split[i] <- TRUE
        info
}
  
me.to.ape.prosse.multi <- function(info, trait.name) {
    if ( nrow(info) == 0 )
      return(NULL)
    Nnode <- sum(!info$split) - 1
    n.tips <- sum(!info$split)

    info$idx2 <- NA
    info$idx2[!info$split] <- 1:n.tips
    info$idx2[as.logical(info$split)] <- order(info$idx[info$split]) + n.tips + 1

    i <- match(info$parent, info$idx)
    info$parent2 <- info$idx2[i]
    info$parent2[is.na(info$parent2)] <- n.tips + 1

    tip.label <- ifelse(subset(info, !split)$extinct,
                        sprintf("ex%d", 1:n.tips),
                        sprintf("sp%d", 1:n.tips))
    node.label <- sprintf("nd%d", 1:Nnode)

    info$name <- NA
    info$name[!info$split] <- tip.label
    
    tmp <-match(1:n.tips, info$idx2)
    
    states <- info$state[tmp]
    names(states) <- tip.label
    
    species <- info$species[tmp]
    names(species) <- tip.label
    
    traits <- as.matrix(info[tmp,trait.name],nrow=n.tips)
    rownames(traits) <- tip.label
    
    phy <- reorder(structure(list(edge=cbind(info$parent2, info$idx2),
                                 Nnode=Nnode,
                                 tip.label=tip.label,
                                 states=states,
                                 species=species,
                                 traits=traits,
                                 node.label=node.label,
                                 edge.length=info$len,
                                 orig=info),
                            class="phylo"))

    phy$edge.state <- info$state[match(phy$edge[,2], info$idx2)]
    phy
}
  
prune.prosse.multi <- function (phy, to.drop = NULL)
{
    if (is.null(to.drop))
        to.drop <- subset(phy$orig, !split)$extinct
    if (sum(!to.drop) < 2) {
        NULL
    }
    else if (any(as.logical(to.drop))) {
        phy2 <- drop.tip.fixed(phy, phy$tip.label[as.logical(to.drop)])
        phy2$orig <- phy2$orig[!phy2$orig$extinct, ]
        phy2$species <- phy2$species[!to.drop]
        phy2$states <- phy2$states[!to.drop]
        phy2$traits <- as.matrix(phy2$traits[!to.drop,])
        phy2$hist <- prune.hist(phy, phy2)
        phy2
    }
    else {
        phy
    }
}

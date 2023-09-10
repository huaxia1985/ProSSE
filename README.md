# ProSSE
Protracted Speciation Extinction

See details in Hua X., Herdha T., Burden C. 2022. Protracted speciation under the state-dependent speciation and extinction approach. Syst. Biol. doi:10.1093/sysbio/syac041

To use ProSSE on lineage-level trees, or slower ode solver for ProSSE on species-level trees, users use the following steps:
1. download all the files in the 'R' folder, except for 'simulation.R'
3. install and library R packages 'diversitree', 'ape', and 'deSolve' in R
4. import all the downloaded R files into R, using source("file location")
5. now all the functions are ready to use

Alternatively, to use faster ode solver for ProSSE on species-level trees, user needs to recompile 'diversitree' R package locally, following these steps:
1. download and upzip 'diversitree' source code
2. download all the files from this repository
3. replace NAMESPACE in the 'diversitree' source code
4. put all the files in the 'R' folder to the 'R' folder in the the 'diversitree' source code
5. put all the files in the 'man' folder to the 'man' folder in the the 'diversitree' source code
6. install R package 'devtools' (note that 'diversitree' needs gsl installed)
7. using devtools to compile 'diversitree' locally

       library(devtools)
       install("your path to save the 'diversitree' source code")
       
8. R may ask if dependent packages need to be updated. Choose according to your needs
9. after the package installed. Restart R.
10. now ProSSE is ready to use 

        library(diversitree)
        
12. ?make.prosse gives details on ProSSE for lineage-level tree
13. ?make.prosse.sp gives details on ProSSE for species-level tree
14. ?tree.prosse gives details on simulating a ProSSE tree


Here is the code you can use to regenerate the results in the paper:

#setting parameter values

    b.list <- seq(0.3,0.6,0.1)
    mu.list <- seq(0,0.2,0.1)
    lambda.list <- c(0.1,0.3,1)
    pars.list <- expand.grid(b.list,mu.list,lambda.list)
    pars.list <- rbind(pars.list,expand.grid(0.7,mu.list[-1],lambda.list))

#setting total simulation time for trees

    max.t=15

#the code simulates 1000 trees under each parameter set. It takes very long time, so users can just simulate 1 tree by changing j in 1:1000 to j in 1.

    n <- nrow(pars.list)
    tree.list <- vector("list",n)
    tree.rep.list <- vector("list",n)
    fit.list <- vector("list",n)
    fit.rep.list <- vector("list",n)
    for (i in 1:n) {
        tree.list[[i]] <- vector("list",1000)
        tree.rep.list[[i]] <- vector("list",1000) 
        fit.list[[i]] <- vector("list",1000)
        fit.rep.list[[i]] <- vector("list",1000)

        for (j in 1:1000) {
        
 #simulate lineage-level tree
 
    pars <- as.numeric(pars.list[i,])
    tree <- try(tree.prosse(pars=pars,max.t=max.t),silent=T)
    while(inherits(tree,"try-error") || is.null(tree) || length(unique(tree$species)) <3 || length(unique(tree$species))>10000) {tree <- try(tree.prosee(pars,max.t),silent=T)}
    tree.list[[i]][[j]] <- tree
            
#get ML estimates using ProSSE for lineage-level tree

    lik <- make.prosse(tree)
    p <- starting.point.prosse(tree)
    fit <- find.mle(lik,p,lower=c(0,0,0))
    fit.list[[i]][[j]] <- fit$par
    
#sample one lineage per species to get species-level tree

    species.label <- unique(tree$species)
    n.species <- length(species.label)
    tip.species <- sapply(species.label,function(i) is.element(tree$species,i))
    n.tip.species <- colSums(tip.species)
    keep.tips <- sapply(1:n.species,function (i) ifelse(n.tip.species[i]==1,which(tip.species[,i]),sample(which(tip.species[,i]),1)))
    to.drop <- rep(TRUE,length(tree$tip.label))
    to.drop[keep.tips] <- FALSE
    tree <- prune.prosse(tree,to.drop)
    tree.rep.list[[i]][[j]] <- tree
            
#get ML estimates using ProSSE for species-level tree

	lik <- make.prosse.sp(tree)
	p <- starting.point.prosse(tree)
	fit <- find.mle(lik,p,lower=c(0,0,0))
	fit.rep.list[[i]][[j]] <- fit$par
	}
 
# Extended ProSSE
Extending Protracted Speciation Extinction to account for multiple modes and rates, and uncertain species identities

See details in Hua X., Moritz C. 2023. A phylogenetic approach to delimitate species in a probabilistic way.

To use the code, the diversitree package needs to be recomplied as described above.

Here is the code you can use to regenerate the simulation results in the paper:

#setting parameter values

	b.mu.list <- rbind(c(0.3,0),c(0.4,0.1),c(0.5,0.1))
	lambda.list <- list(function (x) {0.1*exp(x)}, function (x) {1*exp(x)})
	char <- make.brownian.with.drift(0,0.1)
	trait.q.list <- list(matrix(c(0.1,0.03),nrow=2),matrix(c(0.1,0.1),nrow=2),matrix(c(0.03,0.1),nrow=2))
	trait.p.list <- list(matrix(c(0,0,1,0),nrow=2),matrix(c(0,0.5,0.5,0),nrow=2))
	states.sd <- 0.01

	pars.list <- vector("list",30)
	pars.list[[1]] <- list(b.mu.list[1,1], b.mu.list[1,2], lambda.list[[1]],char,trait.q.list[1],trait.p.list[1])
	pars.list[[2]] <- list(b.mu.list[1,1], b.mu.list[1,2], lambda.list[[1]],char,trait.q.list[1],trait.p.list[2])
	pars.list[[3]] <- list(b.mu.list[1,1], b.mu.list[1,2], lambda.list[[1]],char,trait.q.list[2],trait.p.list[1])
	pars.list[[4]] <- list(b.mu.list[1,1], b.mu.list[1,2], lambda.list[[1]],char,trait.q.list[2],trait.p.list[2])
	pars.list[[5]] <- list(b.mu.list[1,1], b.mu.list[1,2], lambda.list[[1]],char,trait.q.list[3],trait.p.list[1])

	j <- 2
	for (i in 1:5) {
		pars.list[[i+5*(j-1)]] <- pars.list[[i]]
		pars.list[[i+5*(j-1)]][[3]] <- lambda.list[[j]]
	}
	for (j in 2:3) {
	for (i in 1:10) {
		pars.list[[i+10*(j-1)]] <- pars.list[[i]]
		pars.list[[i+10*(j-1)]][[1]] <- b.mu.list[j,1]
		pars.list[[i+10*(j-1)]][[2]] <- b.mu.list[j,2]
	}
	}

#setting total simulation time and size for trees

	max.taxa <- 500
	max.t <- 15

#initiate simuating 100 trees under each parameter set and fit extended ProSSE on each tree using ML

	n <- length(pars.list)
	tree.list <- vector("list",n)
	fit.list <- vector("list",n)
 	mcmc.list <- vector("list",n)

#suppose we use the first parameter set

	i <- 1
	pars <- pars.list[[i]]
	tree.list[[i]] <- vector("list",100)
	fit.list[[i]] <- vector("list",100)
	mcmc.list[[i]] <- vector("list",10)

#get the probability distribution of the root state

	root.p <- as.numeric(pars.list[[i]][[5]][[1]])
	root.p <- rev(root.p/sum(root.p))

#simulate 100 tree under extended ProSSE

	for (j in 1:100) {
 
#sample root state and simulate a tree

	x0 <- sample(c(1,2),1,T,root.p)
	x0 <- c(0,x0)
	tree <- try(tree.prosse.multi(pars,max.taxa,max.t,x0),silent=T)
 
#simulate a new tree if the old tree has too few species, too few lineages, and two few species with each state

	while(inherits(tree,"try-error") || is.null(tree) || length(tree$tip.label)<10 || length(unique(tree$species))<3 || length(unique(tree$species))/length(tree$tip.label)>0.7 || sum(tree$traits==1)/length(tree$tip.label)>0.9 || sum(tree$traits==2)/length(tree$tip.label)>0.9 ) {
        tree <- try(tree.prosse.multi(pars,max.taxa,max.t,x0),silent=T)
	}
	tree.list[[i]][[j]] <- tree

#get ML estimates using extended ProSSE for the tree

	lik <- make.prosse.multi(tree, tree$traits, tree$states, states.sd, exp.x)
	p <- starting.point.prosse.multi(tree, lik, q.div=5, tree$states, states.sd)
 
#you can apply constrains on parameters, e.g., make drift term equals 0

	lik <- constrain(lik,drift~0)
	fit <- find.mle(lik,p[-8],lower=rep(0,8),upper=c(Inf,Inf,Inf,Inf,Inf,1,1,Inf))
	fit.list[[i]][[j]] <- fit$par.full
	}

#pick 10% tips to have unknown species identities

	n.unknown.tip <- round(length(tree$tip.label)*0.1)
	unknown.tip <- sample(x=tree$tip.label,size=n.unknown.tip)
	tree$species[unknown.tip] <- 0
	species <- unique(tree$species)
	tree$species[unknown.tip] <- c(1:n.unknown.tip)+max(species)
	species.name <- tree$species[unknown.tip]

#apply MCMC to infer species identities of these tips, while coestimating parameters

#set prior

	prior <- list(make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.beta(0.5),make.prior.beta(0.5),make.prior.exponential(0.5))
 
#define likelihood function

	lik <- make.prosse.multi(tree, tree$traits, tree$states, states.sd, exp.x)
 	p <- starting.point.prosse.multi(tree, lik, q.div=5, tree$states, states.sd)
	names(p) <- c("b","mu","l.r","q12.1","q21.1","p12.1","p21.1","drift","diffusion")
	x.init <- c(p[-8],tree$species[unknown.tip])
	lik <- constrain(lik,drift~0)
 
#run MCMC

	mcmc.list[[i]][[j]] <- mcmc.prosse.multi(lik=lik, tree=tree, species.name=species.name, unknown.tip=unknown.tip, traits=tree$traits, states=tree$states, states.sd=states.sd, lambda=exp.x, x.init=x.init, nstepsw=30,nsteps=5000,w=rep(1,8),prior=prior,lower=rep(0,8),upper=c(Inf,Inf,Inf,Inf,Inf,1,1,Inf))

#we can also apply MCMC to trees with constant speciation modes and rates that are simulated from ProSSE

	names(tree$species) <- tree$tip.label
	n.unknown.tip <- round(length(tree$tip.label)*0.1)
	unknown.tip <- sample(x=tree$tip.label,size=n.unknown.tip)
	tree$species[unknown.tip] <- 0
	species <- unique(tree$species)
	tree$species[unknown.tip] <- c(1:n.unknown.tip)+max(species)
	species.name <- tree$species[unknown.tip]
	prior <- make.prior.exponential(0.5)
	lik <- make.prosse(tree)
	p <- starting.point.prosse(tree, lik, q.div=5)
	names(p) <- c("b","mu","lambda")
	x.init <- c(p,tree$species[unknown.tip])
	mcmc.list[[i]][[j]] <- mcmc.prosse(lik=lik, tree=tree, species.name=species.name, unknown.tip=unknown.tip, x.init=x.init, nstepsw=30,nsteps=5000,w=rep(1,3),prior=prior,lower=rep(0,3),upper=c(Inf,Inf,Inf))

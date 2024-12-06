# ProSSE
Protracted Speciation Extinction

See details in Hua X., Herdha T., Burden C. 2022. Protracted speciation under the state-dependent speciation and extinction approach. Syst. Biol. doi:10.1093/sysbio/syac041

Extending Protracted Speciation Extinction to account for multiple delimitation criteria, and uncertain species identities

See details in Hua X., Moritz C. 2023. A phylogenetic approach to delimitate species in a probabilistic way.

To use ProSSE without re-compile the diversitree package, users use the following steps (note that this option only allows for solvers in R):
1. download all the files in the 'R' folder
3. install and library R packages 'diversitree', 'ape', and 'deSolve' in R
4. import all the downloaded R files into R, using source("file location")
5. now all the functions are ready to use

Alternatively, to use faster solver in C, user needs to recompile 'diversitree' R package locally, following these steps:
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

Here is the code you can use to regenerate the simulation results in Hua et al. 2022:

setting parameter values

    b.list <- seq(0.3,0.6,0.1)
    mu.list <- seq(0,0.2,0.1)
    lambda.list <- c(0.1,0.3,1)
    pars.list <- expand.grid(b.list,mu.list,lambda.list)
    pars.list <- rbind(pars.list,expand.grid(0.7,mu.list[-1],lambda.list))

setting total simulation time for trees

    max.t=15

the code simulates 1000 trees under each parameter set. It takes very long time, so users can just simulate 1 tree by changing j in 1:1000 to j in 1.

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
    p <- starting.point.prosse(tree,lik)
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
	p <- starting.point.prosse(tree,lik)
	fit <- find.mle(lik,p,lower=c(0,0,0))
	fit.rep.list[[i]][[j]] <- fit$par
	}
 
Here is the code you can use to regenerate the simulation results in Hua and Moritz 2023:

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
#set control=list(method="fftR") if using R solver without recompiling package, which will be slow. Otherwise, the default is using C solver.

	lik <- make.prosse(tree=tree, traits=tree$traits, states=tree$states, states.sd=states.sd, lambda=exp.x, control=list(method="fftR"))
	p <- starting.point.prosse(tree=tree, lik=lik, q.div=5, states=tree$states, states.sd=states.sd, lambda=exp.x)
 
#apply constrains on parameters, e.g., make drift term equals 0

	lik <- constrain(lik,drift~0)
	fit <- find.mle(lik,p[-8],lower=rep(0,8),upper=c(Inf,Inf,Inf,Inf,Inf,1,1,Inf))
	fit.list[[i]][[j]] <- fit$par.full

#apply MCMC to infer species identities of these tips, while coestimating parameters
#pick 10% tips to have unknown species identities

	n.unknown.tip <- round(length(tree$tip.label)*0.1)
	unknown.tip <- sample(x=tree$tip.label,size=n.unknown.tip)

#label species name with numbers

	tree$species[unknown.tip] <- 0
	species <- unique(tree$species)
	tree$species[unknown.tip] <- c(1:n.unknown.tip)+max(species)
	species.name <- tree$species[unknown.tip]

#set prior

	prior <- list(make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.beta(0.5),make.prior.beta(0.5),make.prior.exponential(0.5))
 
#define likelihood function

	x.init <- c(p[-8],tree$species[unknown.tip])
 
#run MCMC
#set control=list(method="fftR") if using R solver without recompiling package, which will be slow. Otherwise, the default is using C solver.

	mcmc.list[[i]][[j]] <- mcmc.prosse(lik=lik, drift~0, tree=tree, species.name=species.name, unknown.tip=unknown.tip, traits=tree$traits, states=tree$states, states.sd=states.sd, lambda=exp.x, control=list(method="fftR"), x.init=x.init, nstepsw=30,nsteps=5000,w=rep(1,8),prior=prior,lower=rep(0,8),upper=c(Inf,Inf,Inf,Inf,Inf,1,1,Inf))
}

#The code you can use to regenerate the case study results in Hua and Moritz 2023 is in the case study folder.
#Here uses the Carlia example to demonstrate how to apply ProSSE to your own data:

#load tree and tip data

    load(Carlia.Rdata)

#the Carlia data includes a sample of trees, let's use the first tree as the example.

    tree <- trees[[1]]

#the tree include some tips that are not Australian Carlia, so we remove them

    tree$tip.label <- gsub("-","_",tree$tip.label)
    droptip <- c(which(substr(tree$tip.label,1,6)!="Carlia"), which(tree$tip.label=="Carlia_bicarinata"),which(tree$tip.label=="Carlia_schlegeli"), which(tree$tip.label=="Carlia_cf_decora"))
    tree <- drop.tip(tree,tip=tree$tip.label[droptip])
    tree$species <- tree$species[-droptip]

#define states, states.sd, and lambda if you want to model speciation completion rate for cryptic species as a function (lambda) of a continuous trait (states). In the Carlia example, the continuous trait is the average temperate over the distribution range of a tip lineage.

    rownames(data) <- data[,1]
    states <- data[tree$species,"state"]
    states.sd <- data[tree$species,"state.sd"]
    lambda <- exp.x
    names(states) <- tree$tip.label
    states <- states-mean(states)
    names(states.sd) <- tree$tip.label
    
#define traits if you want to model occurrence between speciation completion rate and trait transitions. You can include as many traits as you want. In the Carlia example, traits are male throat color and habitat.

    traits <- data[tree$species,c(3,4)]
    rownames(traits) <- tree$tip.label

#define tips with unknown species identities

    unknown.tip <- c("Carlia_vivax_north","Carlia_cf_dogare","Carlia_jarnoldae_cy","Carlia_gracilis_maningrida","Carlia_gracilis_ete","Carlia_gracilis_melville","Carlia_amax_eci","Carlia_amax_gulf","Carlia_amax_ete","Carlia_rufilatus_eci","Carlia_munda_melville","Carlia_storri_cy_ng","Carlia_storri_tv_cy","Carlia_schmeltzii_meq","Carlia_schmeltzii_seq","Carlia_cf_sexdentata_te","Carlia_rhomboidalis_north")
    
#label species name with numbers

    names(tree$tip.label) <- tree$tip.label
    tree$species[unknown.tip] <- tree$tip.label[unknown.tip]
    n.species.max <- length(unique(tree$species))
    species.name <- c(1:n.species.max)
    names(species.name) <- unique(tree$species)
    species <- species.name[tree$species]
    tree$species <- species

#prepare inputs to define likelihood function
#assigning tip with unknown species identities into a species complex if you want to estimate seperate speciation completion rate for cryptic species and morphological species

    types <- numeric(length(tree$tip.label))
    names(types) <- tree$tip.label
    types[c("Carlia_vivax_north","Carlia_vivax_south")] <- 1
    types[c("Carlia_dogare_cy","Carlia_cf_dogare")] <- 2
    types[c("Carlia_jarnoldae_eiu","Carlia_jarnoldae_cy")] <- 3
    types[c("Carlia_gracilis_kim","Carlia_gracilis_wte","Carlia_gracilis_melville","Carlia_gracilis_ete","Carlia_gracilis_maningrida")] <- 4
    types[c("Carlia_amax_ete","Carlia_amax_gulf","Carlia_amax_eci","Carlia_amax_wtek")] <- 5
    types[c("Carlia_isostriacantha","Carlia_triacantha")] <- 6
    types[c("Carlia_insularis","Carlia_johnstonei")] <- 7
    types[c("Carlia_rufilatus_te","Carlia_rufilatus_kim","Carlia_rufilatus_eci")] <- 8
    types[c("Carlia_munda_melville","Carlia_munda_broad","Carlia_munda_ete")] <- 9
    types[c("Carlia_storri_tv_cy","Carlia_storri_cy_ng","Carlia_storri_aru")] <- 10
    types[c("Carlia_schmeltzii_seq","Carlia_schmeltzii_meq","Carlia_schmeltzii_tv_cy")] <- 11
    types[c("Carlia_cf_sexdentata_te","Carlia_sexdentata_cy")] <- 12
    types[c("Carlia_rubrigularis_south","Carlia_crypta")] <- 13
    types[c("Carlia_rhomboidalis_south","Carlia_rhomboidalis_north")] <- 14

#define probable species identities for each tip with unknown species identity, if you have prior information on that. If unknown.list is not provided, then all species identities over the tree plus a unique species name for the tip are used as probable species identities.

    unknown.list <- vector("list",length(unknown.tip)) 
    names(unknown.list) <- unknown.tip
    unknown.list$Carlia_vivax_north <- c("Carlia_vivax","Carlia_vivax_north")
    unknown.list$Carlia_cf_dogare <- c("Carlia_dogare","Carlia_cf_dogare")
    unknown.list$Carlia_jarnoldae_cy <- c("Carlia_jarnoldae","Carlia_jarnoldae_cy")
    unknown.list$Carlia_gracilis_maningrida <- c("Carlia_gracilis","Carlia_gracilis_maningrida","Carlia_gracilis_ete","Carlia_gracilis_melville")
    unknown.list$Carlia_gracilis_ete <- c("Carlia_gracilis","Carlia_gracilis_maningrida","Carlia_gracilis_ete","Carlia_gracilis_melville")
    unknown.list$Carlia_gracilis_melville <- c("Carlia_gracilis","Carlia_gracilis_maningrida","Carlia_gracilis_ete","Carlia_gracilis_melville")
    unknown.list$Carlia_amax_eci <- c("Carlia_amax","Carlia_amax_eci","Carlia_amax_gulf","Carlia_amax_ete")
    unknown.list$Carlia_amax_gulf <- c("Carlia_amax","Carlia_amax_eci","Carlia_amax_gulf","Carlia_amax_ete")
    unknown.list$Carlia_amax_ete <- c("Carlia_amax","Carlia_amax_eci","Carlia_amax_gulf","Carlia_amax_ete")
    unknown.list$Carlia_rufilatus_eci <- c("Carlia_rufilatus_eci","Carlia_rufilatus")
    unknown.list$Carlia_munda_melville <- c("Carlia_munda","Carlia_munda_melville")
    unknown.list$Carlia_storri_cy_ng <- c("Carlia_storri","Carlia_storri_cy_ng","Carlia_storri_tv_cy")
    unknown.list$Carlia_storri_tv_cy <- c("Carlia_storri","Carlia_storri_cy_ng","Carlia_storri_tv_cy")
    unknown.list$Carlia_schmeltzii_meq <- c("Carlia_schmeltzii","Carlia_schmeltzii_meq","Carlia_schmeltzii_seq")
    unknown.list$Carlia_schmeltzii_seq <- c("Carlia_schmeltzii","Carlia_schmeltzii_meq","Carlia_schmeltzii_seq")
    unknown.list$Carlia_cf_sexdentata_te <- c("Carlia_sexdentata","Carlia_cf_sexdentata_te","Carlia_longipes","Carlia_quinquecarinata")
    unknown.list$Carlia_rhomboidalis_north <- c("Carlia_rhomboidalis","Carlia_rhomboidalis_north")

#set priors

    prior <- make.prior.exponential(0.5)

#define likelihood function

    lik <- make.prosse(tree,traits=traits,types=types,states=states,states.sd=states.sd,lambda=lambda)
    lik <- constrain(lik,p12.1~0,p21.1~1,p12.2~1,p21.2~1)
    p <- starting.point.prosse(tree, lik, q.div=5, types=types,lambda=lambda)
    p <- p[1:8]

#set initial parameter values

    x.init <- c(p,species[unknown.tip])
    
#set path to a csv file to save MCMC samples

    filename  <- "~/Carlia result.csv"

#run MCMC

    mcmc.result <- mcmc.prosse(lik=lik, p12.1~0,p21.1~1,p12.2~1,p21.2~1, tree=tree, species.name=species.name, unknown.tip=unknown.tip, unknown.list=unknown.list, traits=traits, types=types, states=states, states.sd=states.sd, lambda=lambda, x.init=x.init, nstepsw=30,nsteps=1000,w=rep(1,8),prior=prior,lower=rep(0,8),upper=rep(Inf,8), save.file=filename)

#define alternative likelihood function, if you want to compare different models. In the Carlia example, the alterantive model is speciation competion rate for cryptic species is constant

    lik.const <- make.prosse(tree,traits=traits,types=types)
    lik.const <- constrain(lik.const,p12.1~0,p21.1~1,p12.2~1,p21.2~1)
    p <- starting.point.prosse(tree, lik.const, q.div=5, types=types)
    p <- p[1:8]
    
#set initial parameter values

    x.init <- c(p,species[unknown.tip])
        
#set path to a csv file to save MCMC samples

    filename  <- "~/Carlia result const.csv"

#run MCMC

    mcmc.result.const <- mcmc.prosse(lik=lik, tree=tree, species.name=species.name, unknown.tip=unknown.tip, unknown.list=unknown.list, traits=traits, types=types, x.init=x.init, nstepsw=30,nsteps=1000,w=rep(1,8),prior=prior,lower=rep(0,8),upper=rep(Inf,8), save.file=filename)



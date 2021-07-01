# ProSSE
Protracted Speciation Extinction

See details in Hua et al. 2021. Protracted speciation under the state-dependent speciation and extinction approach. bioRxiv doi:10.1101/2021.06.29.450466

To implement ProSSE, user needs to recompile 'diversitree' R package locally, following these steps:
1. download and upzip 'diversitree' source code
2. download all the files from this repository
3. replace NAMESPACE in the 'diversitree' source code
4. put all the files in the 'R' folder to the 'R' folder in the the 'diversitree' source code
5. put all the files in the 'man' folder to the 'man' folder in the the 'diversitree' source code
6. install R package 'devtools'
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
    tree <- try(tree.prosse(pars,max.t),silent=T)
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
    }


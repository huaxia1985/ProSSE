rm(list=ls())

#run MCMC on each of 100 trees, each saving in file case_tree number
library(diversitree)
library(deSolve)

load('carlia.Rdata')
set.seed(100)
rownames(data) <- data[,1]
ntrees <- length(trees)
#randomly sample 100 trees
strees <- sample.int(ntrees,100)

for (ii in 1:100) {
mcmc.list <- vector("list",100)
mcmc.const.list <- vector("list",100)
tree <- trees[[strees[ii]]]
#removing tips that are not Australian Carlia
tree$tip.label <- gsub("-","_",tree$tip.label)
droptip <- c(which(substr(tree$tip.label,1,6)!="Carlia"), which(tree$tip.label=="Carlia_bicarinata"),which(tree$tip.label=="Carlia_schlegeli"), which(tree$tip.label=="Carlia_cf_decora"))
tree <- drop.tip(tree,tip=tree$tip.label[droptip])
tree$species <- tree$species[-droptip]

#get trait data
states <- data[tree$species,"state"]
states.sd <- data[tree$species,"state.sd"]
traits <- data[tree$species,c(3,4)]
names(states) <- tree$tip.label
states <- states-mean(states)
names(states.sd) <- tree$tip.label
rownames(traits) <- tree$tip.label
names(tree$species) <- tree$tip.label
names(tree$tip.label) <- tree$tip.label

#define tips with unknown species identities
unknown.tip <- c("Carlia_vivax_north","Carlia_cf_dogare","Carlia_jarnoldae_cy","Carlia_gracilis_maningrida","Carlia_gracilis_ete","Carlia_gracilis_melville","Carlia_amax_eci","Carlia_amax_gulf","Carlia_amax_ete","Carlia_rufilatus_eci","Carlia_munda_melville","Carlia_storri_cy_ng","Carlia_storri_tv_cy","Carlia_schmeltzii_meq","Carlia_schmeltzii_seq","Carlia_cf_sexdentata_te","Carlia_rhomboidalis_north")

#assign each tip to a species complex
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

#label species by numbers
tree$species[unknown.tip] <- tree$tip.label[unknown.tip]
n.species.max <- length(unique(tree$species))
species.name <- c(1:n.species.max)
names(species.name) <- unique(tree$species)
species <- species.name[tree$species]
tree$species <- species

#define probable species identities for each tip
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

#set prior
prior <- make.prior.exponential(0.5)

#define likelihood function
lambda <- exp.x
lik <- make.prosse(tree,traits=traits,types=types,states=states,states.sd=states.sd,lambda=lambda)
lik <- constrain(lik,p12.1~0,p21.1~1,p12.2~1,p21.2~1)

#set initial values
p <- starting.point.prosse(tree, lik, q.div=5, types=types, lambda=lambda)
p <- p[1:8]
x.init <- c(p,species[unknown.tip])

#set filename
filename  <- paste0("~/delimitres",ii,".csv")

#run MCMC
mcmc.list[[ii]] <- 
mcmc.prosse(lik=lik, p12.1~0,p21.1~1,p12.2~1,p21.2~1, tree=tree, species.name=species.name, unknown.tip=unknown.tip, unknown.list=unknown.list, traits=traits, types=types, states=states, states.sd=states.sd, lambda=lambda, x.init=x.init, nstepsw=30,nsteps=1000,w=rep(1,8),prior=prior,lower=rep(0,8),upper=rep(Inf,8), save.file=filename)

#define alternative likelihood function with constant speciation completion rate
lik.const <- make.prosse(tree,traits=traits,types=types)
lik.const <- constrain(lik.const,p12.1~0,p21.1~1,p12.2~1,p21.2~1)

#set initial values
p <- starting.point.prosse(tree, lik.const, q.div=5, types=types)
p <- p[1:8]
x.init <- c(p,species[unknown.tip])

#set filename
filename  <- paste0("~/delimitres.const",ii,".csv")

#run MCMC
mcmc.const.list[[ii]] <-
mcmc.prosse(lik=lik, tree=tree, species.name=species.name, unknown.tip=unknown.tip, unknown.list=unknown.list, traits=traits, types=types, x.init=x.init, nstepsw=30,nsteps=1000,w=rep(1,8),prior=prior,lower=rep(0,8),upper=rep(Inf,8), save.file=filename)
}

#summarise delimitation results
library(diversitree)
library(phytools)
library(deSolve)

theta <- NULL
pp <- NULL
species.list <- NULL

theta.const <- NULL
pp.const <- NULL
species.list.const <- NULL

load('carlia.Rdata')
set.seed(100)
rownames(data) <- data[,1]
ntrees <- length(trees)
strees <- sample.int(ntrees,100)

for (ii in 1:100) {
	tree <- trees[[strees[ii]]]
    #removing tips that are not Australian Carlia
    tree$tip.label <- gsub("-","_",tree$tip.label)
    droptip <- c(which(substr(tree$tip.label,1,6)!="Carlia"), which(tree$tip.label=="Carlia_bicarinata"),which(tree$tip.label=="Carlia_schlegeli"), which(tree$tip.label=="Carlia_cf_decora"))
    tree <- drop.tip(tree,tip=tree$tip.label[droptip])
    tree$species <- tree$species[-droptip]

    #get trait data
    states <- data[tree$species,"state"]
    states.sd <- data[tree$species,"state.sd"]
    traits <- data[tree$species,c(3,4)]
    names(states) <- tree$tip.label
    states <- states-mean(states)
    names(states.sd) <- tree$tip.label
    rownames(traits) <- tree$tip.label
    names(tree$species) <- tree$tip.label
    names(tree$tip.label) <- tree$tip.label

    #define tips with unknown species identities
    unknown.tip <- c("Carlia_vivax_north","Carlia_cf_dogare","Carlia_jarnoldae_cy","Carlia_gracilis_maningrida","Carlia_gracilis_ete","Carlia_gracilis_melville","Carlia_amax_eci","Carlia_amax_gulf","Carlia_amax_ete","Carlia_rufilatus_eci","Carlia_munda_melville","Carlia_storri_cy_ng","Carlia_storri_tv_cy","Carlia_schmeltzii_meq","Carlia_schmeltzii_seq","Carlia_cf_sexdentata_te","Carlia_rhomboidalis_north")

    #assign each tip to a species complex
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

    #label species by numbers
    tree$species[unknown.tip] <- tree$tip.label[unknown.tip]
    n.species.max <- length(unique(tree$species))
    species.name <- c(1:n.species.max)
    names(species.name) <- unique(tree$species)
    species <- species.name[tree$species]
    tree$species <- species
    
    #read MCMC result
    filename <- paste0("~/delimitres",ii,".csv")
    tmp <- read.csv(filename,header=T)
    #discard the first 30 samples
    theta <- rbind(theta,tmp[-c(1:30),2:9])
    pp <- c(pp,tmp[-c(1:30),27])
    species.tmp <- NULL
	for (i in 10:26) {
		species.tmp <- cbind(species.tmp,names(species.name)[tmp[-c(1:30),i]])
	}
	species.list <- rbind(species.list,species.tmp)

    filename <- paste0("~/delimitres.const",ii,".csv")
    tmp <- read.csv(filename,header=T)
    #discard the first 30 samples
    theta.const <- rbind(theta.const,tmp[-c(1:30),2:9])
    pp.const <- c(pp.const,tmp[-c(1:30),27])
    species.tmp <- NULL
    for (i in 10:26) {
        species.tmp <- cbind(species.tmp,names(species.name)[tmp[-c(1:30),i]])
    }
    species.list.const <- rbind(species.list.const,species.tmp)
    }
    colnames(species.list) <- colnames(tmp[,10:26])
    colnames(species.list.const) <- colnames(tmp[,10:26])
}

#count the samples of each possible delimitation for each species complex
table(species.list[,c("Carlia_vivax_north")])
table(species.list[,c("Carlia_cf_dogare")])
table(species.list[,c("Carlia_jarnoldae_cy")])
table(species.list[,c("Carlia_gracilis_melville","Carlia_gracilis_ete","Carlia_gracilis_maningrida")])
table(species.list[,c("Carlia_amax_ete","Carlia_amax_gulf","Carlia_amax_eci")])
table(species.list[,c("Carlia_rufilatus_eci")])
table(species.list[,c("Carlia_munda_melville")])
table(species.list[,c("Carlia_storri_tv_cy","Carlia_storri_cy_ng")])
table(species.list[,c("Carlia_schmeltzii_seq","Carlia_schmeltzii_meq")])
table(species.list[,c("Carlia_cf_sexdentata_te")])
table(species.list[,c("Carlia_rhomboidalis_north")])

#calculating marginal likelihood
library(diversitree)
library(phytools)
library(deSolve)
library(multimode)

##functions needed
getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

getlogVolume <- function (a, eVals) {
    K <- dim(D)[1]
    v <- numeric(K+1)
    v[1] <- log(1)
    v[2] <- log(2)
    for (i in 3:(K+1)) {
        v[i] <- log(2*pi/(i-2))+v[i-2]
    }
    v[K+1]+sum(log(sqrt(eVals)))
}

getEllDraw <- function (theta_star, eVals, eVecs) {
    #theta_star is the vector of mode of the K params
    K <- dim(D)[1]
    rs <- runif(1,min=0,max=1)
    pt <- rnorm(K,mean=0,sd=1) #K*1
    lambda_sur <- sum(pt^2)
    lambda_hyp <- rs^(1/K)/sqrt(lambda_sur)
    ((lambda_hyp*pt)*eVals)%*%t(eVecs) + theta_star #1*K_K*K=1*K
}

geta <- function (theta, pp, theta_star, D) {
    L <- dim(theta)[1]
    K <- dim(theta)[2]
    L_prime <- 0.49*L
    a_min <- a_low <- 0.01
    a_max <- a_high <- 100
    delta_a <- a_high-a_low
    c <- 0
    rho <- median(pp)
    while (delta_a > 0.1) {
        delta_a <- a_high-a_low
        D_inv <- solve(D)
        L_low <- L_high <- 0
        Z <- as.matrix(theta-t(theta_star%*%t(rep(1,length(pp)))))
        L_low <- sum( (pp > rho) * (rowSums((Z%*%((1/a_low)*D_inv))*Z) < 1))
        L_high <- sum( (pp > rho) * (rowSums((Z%*%((1/a_high)*D_inv))*Z) < 1))
        if (L_high >= L_prime) {
            if (c==0) {
                a_max <- max(a_high,a_max)
            } else {
                if (a_high < a_max) {
                    a_max <- a_high
                }
            }
            a_high <- (a_low+a_high)/2
            c <- 1
        } else {
            if (c==1) {
                a_high <- min(a_high + 10, a_max)
            } else {
                a_high <- a_high + 100
            }
        }
        if (L_low < L_prime) {
            if (a_low > a_min) {
                a_min <- a_low
            }
            a_low <- (a_low + a_high)/2
        } else {
            a_low <- max((1+a_low)/2,a_min)
        }
    }
    a_high
}

#marginal likeihoods from the posterior draws through a geometric identity
#gives log marginal likelihood
#theta is posterior draws of params
#pp is the log posterior probabiltiy of theta
#R is the large number of draws from A2
#posterior is the function to calculate log posterior probability

calMarLik <- function (theta, pp, R, posterior) {
    rho <- median(pp)
    L <- length(pp)
    theta_star <- apply(theta, 2, getmode)
    D <- cov(theta)
    a <- geta(theta, pp, theta_star, D)
    eigen_aD <- eigen(a*D)
    eVals <- eigen_aD$values
    eVecs <- eigen_aD$vectors
    pp_l <- NULL
    for (i in 1:L) {
            if ((pp[i] > rho) && (((theta[i,]-theta_star)%*%((1/a)*D_inv)%*%t(theta[i,]-theta_star)) < 1)) {
                pp_l <- c(pp_l,pp[i])
            }
    }
    ln_k_A <- -rho+log(sum(exp(-pp_l+rho)))-log(L)
    theta_r <- replicate(R, getEllDraw(theta_star, eVals, eVecs))
    pp_r <- apply(theta_r, 1, posterior)
    r_bar <- sum(pp_r>rho)
    pi_hat <- r-bar/R
    V_A <- pi_hat*getVolumn(a,eVals)
    log(V_A)-log(k_A)
}

rho <- median(pp)
L <- length(pp)
theta_star <- sapply(apply(theta,2,locmodes,1),function (i) i$locations)
D <- cov(theta)
a <- geta(theta, pp, theta_star, D)
eigen_aD <- eigen(a*D)
eVals <- eigen_aD$values
eVecs <- eigen_aD$vectors
D_inv <- solve(D)

Z <- as.matrix(theta-t(theta_star%*%t(rep(1,length(pp)))))
tmp <- which((pp > rho) * (rowSums((Z%*%((1/a)*D_inv))*Z) < 1)==1)
pp_l <- pp[tmp]

ln_k_A <- -rho+log(sum(exp(-pp_l+rho)))-log(L)

r_bar <- 0

load('carlia.Rdata')
set.seed(100)
rownames(data) <- data[,1]
ntrees <- length(trees)
#randomly sample 100 trees
strees <- sample.int(ntrees,100)

for (ii in 1:100) {
	tree <- trees[[strees[ii]]]
    #removing tips that are not Australian Carlia
    tree$tip.label <- gsub("-","_",tree$tip.label)
    droptip <- c(which(substr(tree$tip.label,1,6)!="Carlia"), which(tree$tip.label=="Carlia_bicarinata"),which(tree$tip.label=="Carlia_schlegeli"), which(tree$tip.label=="Carlia_cf_decora"))
    tree <- drop.tip(tree,tip=tree$tip.label[droptip])
    tree$species <- tree$species[-droptip]

    #get trait data
    states <- data[tree$species,"state"]
    states.sd <- data[tree$species,"state.sd"]
    traits <- data[tree$species,c(3,4)]
    names(states) <- tree$tip.label
    states <- states-mean(states)
    names(states.sd) <- tree$tip.label
    rownames(traits) <- tree$tip.label
    names(tree$species) <- tree$tip.label
    names(tree$tip.label) <- tree$tip.label

    #define tips with unknown species identities
    unknown.tip <- c("Carlia_vivax_north","Carlia_cf_dogare","Carlia_jarnoldae_cy","Carlia_gracilis_maningrida","Carlia_gracilis_ete","Carlia_gracilis_melville","Carlia_amax_eci","Carlia_amax_gulf","Carlia_amax_ete","Carlia_rufilatus_eci","Carlia_munda_melville","Carlia_storri_cy_ng","Carlia_storri_tv_cy","Carlia_schmeltzii_meq","Carlia_schmeltzii_seq","Carlia_cf_sexdentata_te","Carlia_rhomboidalis_north")

    #assign each tip to a species complex
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

    #label species by numbers
    tree$species[unknown.tip] <- tree$tip.label[unknown.tip]
    n.species.max <- length(unique(tree$species))
    species.name <- c(1:n.species.max)
    names(species.name) <- unique(tree$species)
    species <- species.name[tree$species]
    tree$species <- species
    
    #set prior
    prior <- make.prior.exponential(0.5)

    #define likelihood function
    lambda <- exp.x
    lik <- make.prosse(tree,traits=traits,types=types,states=states,states.sd=states.sd,lambda=lambda)
    lik <- constrain(lik,p12.1~0,p21.1~1,p12.2~1,p21.2~1)
    
    #define posterior
	posterior <- function (x,prior,tree) {
		if (any(x[1:8]<0)) {
			out <- -Inf
		} else {
			tree$species[unknown.tip] <- x[9:25]
			if (check.paraphyletic(tree)) {
				out <- -Inf
			} else {
				lik <- make.prosse(tree,traits=traits,types=types,states=states,states.sd=states.sd,lambda=lambda)
				lik <- constrain(lik,p12.1~0,p21.1~1,p12.2~1,p21.2~1)
				out <- lik(x[1:8]) + sum(prior(x[1:8]))
        		if(is.na(out) || is.infinite(out)) {
        			out <- -Inf
        		}
        	}
   		}
        out
	}
    
    #calculate posterior of 500 samples from the manifold
	theta_r <- replicate(500, getEllDraw(theta_star, eVals, eVecs),simplify="matrix")
	species_r <- matrix(species.name[species.list[sample.int(n=1000,size=500)+(ii-1)*1000,]],nrow=500)
	x.init_r <- cbind(t(theta_r),species_r)
	colnames(x.init_r) <- c(argnames(lik),species[unknown.tip])

	pp_r <- apply(x.init_r, 1, posterior, prior, tree)
	r_bar <- r_bar + sum(pp_r>rho)
}
pi_hat <- r_bar/500/100
ln_V_A <- log(pi_hat)+getlogVolume(a,eVals)
marginal_lik <- ln_V_A-ln_k_A

#calculating marginal likelihood for constant speciation completion rate
rho <- median(pp.const)
L <- length(pp.const)
theta_star <- sapply(apply(theta.const,2,locmodes,1),function (i) i$locations)
D <- cov(theta.const)
a <- geta(theta.const, pp.const, theta_star, D)
eigen_aD <- eigen(a*D)
eVals <- eigen_aD$values
eVecs <- eigen_aD$vectors
D_inv <- solve(D)

Z <- as.matrix(theta.const-t(theta_star%*%t(rep(1,length(pp.const)))))
tmp <- which((pp.const > rho) * (rowSums((Z%*%((1/a)*D_inv))*Z) < 1)==1)
pp_l <- pp.const[tmp]

ln_k_A <- -rho+log(sum(exp(-pp_l+rho)))-log(L)

r_bar <- 0

for (ii in 1:100) {
    tree <- trees[[strees[ii]]]
    #removing tips that are not Australian Carlia
    tree$tip.label <- gsub("-","_",tree$tip.label)
    droptip <- c(which(substr(tree$tip.label,1,6)!="Carlia"), which(tree$tip.label=="Carlia_bicarinata"),which(tree$tip.label=="Carlia_schlegeli"), which(tree$tip.label=="Carlia_cf_decora"))
    tree <- drop.tip(tree,tip=tree$tip.label[droptip])
    tree$species <- tree$species[-droptip]

    #get trait data
    states <- data[tree$species,"state"]
    states.sd <- data[tree$species,"state.sd"]
    traits <- data[tree$species,c(3,4)]
    names(states) <- tree$tip.label
    states <- states-mean(states)
    names(states.sd) <- tree$tip.label
    rownames(traits) <- tree$tip.label
    names(tree$species) <- tree$tip.label
    names(tree$tip.label) <- tree$tip.label

    #define tips with unknown species identities
    unknown.tip <- c("Carlia_vivax_north","Carlia_cf_dogare","Carlia_jarnoldae_cy","Carlia_gracilis_maningrida","Carlia_gracilis_ete","Carlia_gracilis_melville","Carlia_amax_eci","Carlia_amax_gulf","Carlia_amax_ete","Carlia_rufilatus_eci","Carlia_munda_melville","Carlia_storri_cy_ng","Carlia_storri_tv_cy","Carlia_schmeltzii_meq","Carlia_schmeltzii_seq","Carlia_cf_sexdentata_te","Carlia_rhomboidalis_north")

    #assign each tip to a species complex
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

    #label species by numbers
    tree$species[unknown.tip] <- tree$tip.label[unknown.tip]
    n.species.max <- length(unique(tree$species))
    species.name <- c(1:n.species.max)
    names(species.name) <- unique(tree$species)
    species <- species.name[tree$species]
    tree$species <- species
    
    #set prior
    prior <- make.prior.exponential(0.5)

    #define likelihood function
    lik <- make.prosse(tree,traits=traits,types=types)
    lik <- constrain(lik.const,p12.1~0,p21.1~1,p12.2~1,p21.2~1)
    
    #define posterior
    posterior <- function (x,prior,tree) {
        if (any(x[1:8]<0)) {
            out <- -Inf
        } else {
            tree$species[unknown.tip] <- x[9:25]
            if (check.paraphyletic(tree)) {
                out <- -Inf
            } else {
                lik <- make.prosse(tree,traits=traits,types=types)
                lik <- constrain(lik,p12.1~0,p21.1~1,p12.2~1,p21.2~1)
                out <- lik(x[1:8]) + sum(prior(x[1:8]))
                if(is.na(out) || is.infinite(out)) {
                    out <- -Inf
                }
            }
           }
        out
    }
    
    #calculate posterior of 500 samples from the manifold
    theta_r <- replicate(500, getEllDraw(theta_star, eVals, eVecs),simplify="matrix")
    species_r <- matrix(species.name[species.list.const[sample.int(n=1000,size=500)+(ii-1)*1000,]],nrow=500)
    x.init_r <- cbind(t(theta_r),species_r)
    colnames(x.init_r) <- c(argnames(lik),species[unknown.tip])

    pp_r <- apply(x.init_r, 1, posterior, prior, tree)
    r_bar <- r_bar + sum(pp_r>rho)
}
pi_hat <- r_bar/500/100
ln_V_A <- log(pi_hat)+getlogVolume(a,eVals)
marginal_lik_const <- ln_V_A-ln_k_A


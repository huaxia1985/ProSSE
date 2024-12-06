rm(list=ls())

#download lionepha.zip from https://jeetsukumaran.github.io/delineate/workflow1.html
#setwd(where the file is unzipped)

library(diversitree)
library(phytools)
tree <- read.nexus("lionepha/04-species-delimitation/lionepha-p095-hkyg.mcct-mean-age.tree.nex")
run <- read.table("lionepha/04-species-delimitation/lionepha.run1.tsv",sep="\t",header=T)
tree$species <- as.character(run$species)

#run$status=0 indicates lineages with uncertain species identities
tree$species[run$status==0] <- as.character(run$lineage[run$status==0])
names(tree$species) <- run$lineage
tree$species <- tree$species[tree$tip.label]
unknown.tip <- as.character(run$lineage[run$status==0])

#label species with numbers
n.species.max <- length(unique(tree$species))
species.name <- c(1:n.species.max)
names(species.name) <- unique(tree$species)
species <- species.name[tree$species]
tree$species <- species

#set prior
prior <- make.prior.exponential(0.5)

#scale tree length to make rate parameter values not too large
tree$edge.length <- tree$edge.length*10^4
tree$root.depth <- max(nodeHeights(tree))

#define likelihood function
lik <- make.prosse(tree)

#set initial values
p <- starting.point.prosse(tree, lik, q.div=5)
x.init <- c(p,species[unknown.tip])

mcmc.list <- mcmc.prosse(lik=lik, tree=tree, species.name=species.name, unknown.tip=unknown.tip, unknown.list=NULL, x.init=x.init, nstepsw=30,nsteps=5000,w=rep(1,3),prior=prior,lower=rep(0,3),upper=c(Inf,Inf,Inf))
}
#summarise delimitation results, burin 30
table(mcmc.list[-c(1:30),c("coalescentpop008","coalescentpop009")])
#2,3=new species; 1=casta; 7=lindrothi
table(mcmc.list[-c(1:30),c("L_lindrothi_CA_Emerald_Lake","coalescentpop011","L_lindrothi_CA_South_Fork_Bishop_Creek")])
#4,5,6=new species; 7= lindrothi
table(mcmc.list[-c(1:30),c("L_probata_UT_Shingle_Creek","coalescentpop024")]
+ )
#10=probata; 11,12=new species
table(mcmc.list[-c(1:30),c("L_disjuncta_OR_Mt_Hood","L_disjuncta_CA_Lily_Lake","coalescentpop006")]            )
#13,14,15=new species, 16=disjuncta
table(mcmc.list[-c(1:30),c("L_tuulukwa_CA_Trinity_Alps","L_tuulukwa_OR_Knowles_Creek","L_tuulukwa_OR_Marys_Peak")]
+ )
#17=osculans,18,19,20=new species,21=sequoiae

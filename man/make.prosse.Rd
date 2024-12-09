\name{make.prosse}
\alias{make.prosse}
\alias{starting.point.prosse}

\title{Protracted Speciation and Extinction Model for lineage-level tree}

\description{Prepare to run ProSSE (Protracted Speciation and
  Extinction) on a lineage-level tree and species identity of each tip.  This
  function creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.}

\usage{
make.prosse(tree, traits, types, states, states.sd, lambda, control)
starting.point.prosse(tree, lik, q.div=5, types, states, states.sd, lambda)
}

\arguments{
  \item{tree}{An ultrametric bifurcating phylogenetic tree, in
    \code{ape} \dQuote{phylo} format. In the tree, tree$species, 
    named by tree$tip.label, gives the species identity of each tip}
  \item{traits}{A dataframe with each column corresponding to a trait 
  	and each row a tip lineage. Each row is named by tree$tip.label.
  	When provided, speciation completion is allowed to co-occur with
  	trait transition}
  \item{types}{A vector of integers that label the species complex
  	a tip lineage is from. Each element is named by tree$tip.label.
  	When provided, seperate speciation completion rates are estimated
  	for cryptic species and morphological species}
  \item{states}{A vector of the value of a continuous trait that a 
  	tip lineage has. Each element is named by 
  	tree$tip.label. When provided, cryptic speciation completion rate 
  	is modelled as a function of the continuous trait}
  \item{states.sd}{A vector of the standard deviation of the continuous 
  	trait that a tip lineage has. Each element is named by tree$tip.label.
  	When provided, variation in the values of the continuous trait is
  	accounted for}
  \item{lambda}{A function that defines how cryptic speciation completion 
  	rate depends on the continuous trait. Lambda is required when states
  	is provided}
  \item{control}{A list of control parameters for the FFT. See make.quasse
    for details. Need to set control=list(method="fftR") to use solvers in 
    R when using prosse without recompling}
  \item{q.div}{see starting.point.bisse}
}

\details{
  \code{make.prosse} returns a function of class \code{prosse}.

  \preformatted{
    f(pars, condition.surv=TRUE)
  }

  The arguments are interpreted as
  \itemize{
    \item \code{pars} A vector of parameters. When only tree is provided,
    the parameters are \code{b}, \code{mu}, \code{lambda}, which are speciation 
    initiation rate, extinction rate, speciation completion rate. When traits 
    is provided, additional parameters are \code{q} and \code{p} for each possible
    trait transition, which are the transition rate and the probability that 
    speciation completion co-occurs with trait transition. When types is provided,
    parameter \code{lambda} is seperately estimated for cryptic species
    \code{lambda.c} and for morphological species \code{lambda.m}. When states 
    is provided, additional parameters are the coefficients in lambda function and
    drift and diffusion for the continuous trait.
    
    \item \code{condition.surv} (logical): should the likelihood
    calculation condition on survival of two lineages and the speciation
    event subtending them?  This is done by default, following Nee et
    al. 1994.
    }

  \code{starting.point.bisse} produces a heuristic starting point to
  start from.
}

\examples{
#assume a constant speciation completion rate across the tree	
pars <- c(0.3,0.1,0.1)
phy <- tree.prosse(pars,max.t=15)
lik <- make.prosse(phy)
p <- starting.point.prosse(phy,lik)
fit <- find.mle(lik,p,lower=c(0,0,0))

#assume different speciation completion rates for cryptic species and morpological species
load(Carlia.Rdata)
tree <- trees[[1]]
tree$tip.label <- gsub("-","_",tree$tip.label)
droptip <- c(which(substr(tree$tip.label,1,6)!="Carlia"), which(tree$tip.label=="Carlia_bicarinata"),which(tree$tip.label=="Carlia_schlegeli"), which(tree$tip.label=="Carlia_cf_decora"))
tree <- drop.tip(tree,tip=tree$tip.label[droptip])
tree$species <- tree$species[-droptip]
unknown.tip <- c("Carlia_vivax_north","Carlia_cf_dogare","Carlia_jarnoldae_cy","Carlia_gracilis_maningrida","Carlia_gracilis_ete","Carlia_gracilis_melville","Carlia_amax_eci","Carlia_amax_gulf","Carlia_amax_ete","Carlia_rufilatus_eci","Carlia_munda_melville","Carlia_storri_cy_ng","Carlia_storri_tv_cy","Carlia_schmeltzii_meq","Carlia_schmeltzii_seq","Carlia_cf_sexdentata_te","Carlia_rhomboidalis_north")
names(tree$tip.label) <- tree$tip.label
tree$species[unknown.tip] <- tree$tip.label[unknown.tip]
n.species.max <- length(unique(tree$species))
species.name <- c(1:n.species.max)
names(species.name) <- unique(tree$species)
species <- species.name[tree$species]
tree$species <- species
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
lik <- make.prosse(tree,types=types)
p <- starting.point.prosse(tree, lik, types=types)
}

\references{
  Hua X., Herdha T., and Burden C. 2021. Protracted speciation under the
  state-dependent speciation and extinction approach. bioRxiv doi:10.1101/2021.06.29.450466
  Hua X., Moritz C. 2023. A phylogenetic approach to delimitate species in a probabilistic way.
}

\author{Xia Hua}
\keyword{models}

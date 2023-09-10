\name{mcmc.prosse}
\alias{mcmc.prosse}

\title{Infer species indentity using ProSSE}

\description{Coestimate ProSSE parameters using slice sampling and
	species indentities using Gibbs sampling}

\usage{
mcmc.prosse(lik, tree, species.name, unknwon.tip, x.init, nsteps, w, ...)
\method{mcmc.prosse}{default}(lik, tree, species.name, unknwon.tip, unknown.list=NULL, x.init, nstepsw=0, 
			nsteps, w, prior=NULL, sampler=sampler.slice, fail.value=-Inf, lower=-Inf,
             upper=Inf, print.every=1, control=list(),
             save.file, save.every=0, save.every.dt=NULL,
             previous=NULL, previous.tol=1e-4, keep.func=TRUE, ...)
}

\arguments{
  \item{lik}{Likelihood function generated from make.prosse.}
  
  \item{tree}{An ultrametric bifurcating phylogenetic tree, in
    \code{ape} \dQuote{phylo} format. In the tree, tree$species, 
    named by tree$tip.label, gives the species identity of each tip}
    
  \item{species.name}{species name of the tips with unknown species 
  	identities if the tip is an unique species}
  	
  \item{unknown.tip}{tip name of the tips with unknown species identities}
  
  \item{unknown.list}{a list named by unknwon.tip, gives the set of potential
  	species identities of the tip. When there is no prior knowledge on the set,
  	unknown.list is set to NULL}
  	
  \item{x.init}{Initial parameter values.}
  
  \item{nsteps}{Number of MCMC steps to take.}
  
  \item{w}{Tuning parameter for the sampler.  See ?mcmc for details.}

  \item{nstepsw}{number of the initial steps to determine w}
  
}

\details{
	see ?mcmc for details on the other default inputs.
}

\examples{
pars <- c(b=0.4,mu=0.1,lambda=0.1)
tree <- tree.prosse(pars=pars,max.t=15)
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
res <- mcmc.prosse(lik=lik, tree=tree, species.name=species.name, unknown.tip=unknown.tip, x.init=x.init, nstepsw=30,nsteps=5000,w=rep(1,3),prior=prior,lower=rep(0,3),upper=c(Inf,Inf,Inf))
}

\references{
  Hua X., Moritz C. 2023. A phylogenetic approach to delimitate species in a
  probabilistic way.
}

\author{Xia Hua}
\keyword{models}
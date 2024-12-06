\name{mcmc.prosse}
\alias{mcmc.prosse}

\title{Infer species indentity using ProSSE}

\description{Coestimate ProSSE parameters using slice sampling and
	species indentities using Gibbs sampling}

\usage{
mcmc.prosse(lik, ..., tree, species.name, unknwon.tip, unknown.list, 
	x.init, nstepsw, nsteps, w, prior, lower, upper, save.file)
\method{mcmc.prosse}{default}(lik, ..., tree, species.name, unknown.tip, 
	unknown.list=NULL, traits=NULL, types=NULL, states=NULL, states.sd=NULL, 
	lambda=NULL, control=NULL, x.init, nstepsw=0, nsteps, w, prior=NULL,
    lower=-Inf, upper=Inf, sampler=sampler.slice, fail.value=-Inf,
    print.every=1,save.file, save.every=0, save.every.dt=NULL,
    previous=NULL,previous.tol=1e-4,keep.func=TRUE)
}

\arguments{
  \item{lik}{Likelihood function generated from make.prosse. See ?make.prosse
  	for arguments traits, types, states, states.sd, lambda, control}
  
  \item{...}{other arguments, for example, constraints on the 
  	likelihood function}
  
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
  
  \item{nstepsw}{number of the initial steps to determine w}
  
  \item{nsteps}{Number of MCMC steps to take.}
  
  \item{w}{Tuning parameter for the sampler.  See ?mcmc for details.}
  
  \item{prior}{A function if the same prior is used for all parameters.
  	Otherwise, a list of functions, each gives the prior for a corresponding parameter}
  	
  \item{lower}{A vector giving the lower bounds of parameters}
  
  \item{upper}{A vector giving the upper bounds of parameters}
  
  \item{save.file}{Path to a csv file to save each posterior sample}
}

\details{
	see ?mcmc for details on the other default inputs.
}

\examples{
#assuming constant speciation completion rate across the tree
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
x.init <- c(p,tree$species[unknown.tip])
res <- mcmc.prosse(lik=lik, tree=tree, species.name=species.name, unknown.tip=unknown.tip, x.init=x.init, nstepsw=30,nsteps=5000,w=rep(1,3),prior=prior,lower=rep(0,3),upper=c(Inf,Inf,Inf))
#assuming speciation completion rate co-occur with a binary trait and depends on a continuous trait
lambda <- function (x) {0.1*exp(x)}
char <- make.brownian.with.drift(0,0.1)
trait.q <- matrix(c(0.1,0.03),nrow=2)
trait.p <- matrix(c(0,0,1,0),nrow=2)
states.sd <- 1
root.p <- rev(trait.q/sum(trait.q))
x0 <- sample(c(1,2),1,T,root.p)
x0 <- c(0,x0)
pars <- list(0.4,0.1,lambda,char,list(trait.q),list(trait.p))
tree <- tree.prosse.multi(pars,max.taxa=500,max.t=15,x0)
lik <- make.prosse(tree,traits=tree$traits,states=tree$states,states.sd=states.sd,lambda=exp.x,control=list(method="fftR"))
lik <- constrain(lik,drift~0)
p <- starting.point.prosse(tree,lik,states=tree$states,states.sd=states.sd,lambda=exp.x)
names(tree$species) <- tree$tip.label
n.unknown.tip <- round(length(tree$tip.label)*0.1)
unknown.tip <- sample(x=tree$tip.label,size=n.unknown.tip)
tree$species[unknown.tip] <- 0
species <- unique(tree$species)
tree$species[unknown.tip] <- c(1:n.unknown.tip)+max(species)
species.name <- tree$species[unknown.tip]
x.init <- c(p[-8],tree$species[unknown.tip])
prior <- list(make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.beta(0.5),make.prior.beta(0.5),make.prior.exponential(0.5))
res <- mcmc.prosse(lik=lik, drift~0, tree=tree, species.name=species.name, unknown.tip=unknown.tip, traits=tree$traits, states=tree$states, states.sd=states.sd, lambda=exp.x, x.init=x.init, nstepsw=30,nsteps=5000,w=rep(1,8),prior=prior,lower=rep(0,8),upper=c(Inf,Inf,Inf,Inf,Inf,1,1,Inf))
}

\references{
  Hua X., Moritz C. 2023. A phylogenetic approach to delimitate species in a
  probabilistic way.
}

\author{Xia Hua}
\keyword{models}
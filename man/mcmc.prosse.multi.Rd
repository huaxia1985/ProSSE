\name{mcmc.prosse.multi}
\alias{mcmc.prosse.multi}

\title{Infer species indentity using extended ProSSE for multiple speciation modes and rates}

\description{Coestimate extended ProSSE parameters using slice sampling and
	species indentities using Gibbs sampling}

\usage{
mcmc.prosse.multi(lik, tree, species.name, unknwon.tip, unknown.list, traits, states, 
	states.sd, lambda, constraint, x.init, nsteps, w, ...)
\method{mcmc.prosse}{default}(lik, tree, species.name, unknwon.tip, unknown.list, 
	traits, states, states.sd, lambda, constraint, control=NULL, x.init, nstepsw=0, 
			nsteps, w, prior=NULL, sampler=sampler.slice2, fail.value=-Inf, lower=-Inf,
             upper=Inf, print.every=1, save.file, save.every=0, save.every.dt=NULL,
             previous=NULL, ...)
}

\arguments{
  \item{lik}{Likelihood function generated from make.prosse.multi.}
  
  \item{tree}{An ultrametric bifurcating phylogenetic tree, in
    \code{ape} \dQuote{phylo} format. In the tree, tree$species, 
    named by tree$tip.label, gives the species identity of each tip}
    
  \item{species.name}{species name of the tips with unknown species 
  	identities if the tip is an unique species}
  	
  \item{unknown.tip}{tip name of the tips with unknown species identities}
  
  \item{unknown.list}{a list named by unknwon.tip, gives the set of potential
  	species identities of the tip. When there is no prior knowledge on the set,
  	unknown.list is set to NULL}
  
  \item{traits}{A matrix with each row named by tree$tip.label and each
  	column being a discrete trait}
  	
  \item{states}{A vector named by tree$tip.label, gives the value of
  	a continuous trait of each tip; the continuous trait is assumed to 
  	affect speciation completion rate}
  	
  \item{states.sd}{A value or a vector named by tree$tip.label, gives
  	the standard deviation of the continuous trait of each tip}
  	
  \item{lambda}{A function describes how lambda depends on the continuous
  	trait}
  
  \item{constraint}{constraints on parameters}
  	
  \item{control}{A list of parameters for tuning the performance of the
  	integrator. See ?make.quasse for details}

  \item{x.init}{Initial parameter values.}
  
  \item{nsteps}{Number of MCMC steps to take.}
  
  \item{w}{Tuning parameter for the sampler.  See ?mcmc for details.}

  \item{nstepsw}{number of the initial steps to determine w}
  
}

\details{
	see ?mcmc for details on the other default inputs.
}

\examples{
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
lik <- make.prosse.multi2(tree,tree$traits,tree$states,states.sd,exp.x)
lik <- constrain(lik,drift~0)
p <- starting.point.prosse.multi(tree,lik,q.div=5,tree$states,states.sd)
names(p) <- c("b","mu","l.r","q12.1","q21.1","p12.1","p21.1","drift","diffusion")
names(tree$species) <- tree$tip.label
n.unknown.tip <- round(length(tree$tip.label)*0.1)
unknown.tip <- sample(x=tree$tip.label,size=n.unknown.tip)
tree$species[unknown.tip] <- 0
species <- unique(tree$species)
tree$species[unknown.tip] <- c(1:n.unknown.tip)+max(species)
species.name <- tree$species[unknown.tip]
x.init <- c(p[-8],tree$species[unknown.tip])
prior <- list(make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.exponential(0.5),make.prior.beta(0.5),make.prior.beta(0.5),make.prior.exponential(0.5))
res <- mcmc.prosse.multi(lik=lik, tree=tree, species.name=species.name, unknown.tip=unknown.tip, traits=tree$traits, states=tree$states, states.sd=states.sd, lambda=exp.x, constraint=list(drift~0),x.init=x.init, nstepsw=30,nsteps=5000,w=rep(1,8),prior=prior,lower=rep(0,8),upper=c(Inf,Inf,Inf,Inf,Inf,1,1,Inf))
}

\references{
  Hua X., Moritz C. 2023. A phylogenetic approach to delimitate species in a
  probabilistic way.
}

\author{Xia Hua}
\keyword{models}
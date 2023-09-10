\name{make.prosse.multi}
\alias{make.prosse.multi}
\alias{starting.point.prosse}

\title{Extended Protracted Speciation and Extinction Model for multiple modes and rates}

\description{Prepare to run Extended ProSSE on a lineage-level 
  tree and species identity of each tip.  This function creates 
  a likelihood function that can be used in \link[=find.mle]
  {maximum likelihood} or \link[=mcmc]{Bayesian} inference.}

\usage{
make.prosse.multi(tree, traits, states, states.sd, lambda, control=NULL)
make.prosse.multi2(tree, traits, states, states.sd, lambda, control=NULL)
starting.point.prosse.multi(tree, lik, q.div=5, states, states.sd)
}

\arguments{
  \item{tree}{An ultrametric bifurcating phylogenetic tree, in
    \code{ape} \dQuote{phylo} format. In the tree, tree$species, 
    named by tree$tip.label, gives the species identity of each tip}
  \item{traits}{A matrix with each row named by tree$tip.label and each
  	column being a discrete trait}
  \item{states}{A vector named by tree$tip.label, gives the value of
  	a continuous trait of each tip; the continuous trait is assumed to 
  	affect speciation completion rate}
  \item{states.sd}{A value or a vector named by tree$tip.label, gives
  	the standard deviation of the continuous trait of each tip}
  \item{lambda}{A function describes how lambda depends on the continuous
  	trait}
  \item{control}{A list of parameters for tuning the performance of the
  	integrator. See ?make.quasse for details}
}

\details{
  \code{make.prosse.multi} returns a function of class \code{prosse.multi}.
  \code{make.prosse.multi} returns a function of class \code{prosse.multi} 
  with output compatible to mcmc.prosse.multi.
  
  \preformatted{
    f(pars, condition.surv=TRUE, root=ROOT.OBS)
  }

  The arguments are interpreted as
  \itemize{
    \item \code{pars} A vector of parameters, in the order
    of speciation initiation rate, extinction rate, coefficients in the 
    lambda function, transition rates in discrete traits, probability of 
    speciation completion event co-occuring with trait transitions, drift
    and diffusion in the continuous trait.
    \item \code{condition.surv} (logical): should the likelihood
    calculation condition on survival of two lineages and the speciation
    event subtending them?  This is done by default, following Nee et
    al. 1994.
    \item \code{root} Behaviour at the root. See ?make.bisse for details.
    }

  \code{starting.point.prosse.multi} produces a heuristic starting point to
  start from.
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
phy <- tree.prosse.multi(pars,max.taxa=500,max.t=15,x0)
lik <- make.prosse.multi(phy,phy$traits,phy$states,states.sd,exp.x)
p <- starting.point.prosse.multi(phy,lik,q.div=5,phy$states,states.sd)
lik <- constrain(lik,drift~0)
fit <- find.mle(lik,p[-8],lower=rep(0,8),upper=c(Inf,Inf,Inf,Inf,Inf,1,1,Inf))
}

\references{
  Hua X., Moritz C. 2023. A phylogenetic approach to delimitate species in a
  probabilistic way.
}

\author{Xia Hua}
\keyword{models}

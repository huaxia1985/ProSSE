\name{make.prosse.sp}
\alias{make.prosse.sp}

\title{Protracted Speciation and Extinction Model for species-level tree}

\description{Prepare to run ProSSE.sp (Protracted Speciation and
  Extinction for species-level tree) on a tree that includes only 
  one representative lineage of each extant species.  
  This function creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.}

\usage{
make.prosse.sp(tree, control=list())
}

\arguments{
  \item{tree}{An ultrametric bifurcating phylogenetic tree, in
    \code{ape} \dQuote{phylo} format. The tree includes only one representative
  lineage of each extant species. }
  
  \item{control}{List of control parameters for the ODE solver.  See make.bisse
    for details.}
}

\details{
  \code{make.prosse.sp} returns a function of class \code{prosse.sp}.  This
  function has argument list (and default values)

  \preformatted{
    f(pars, condition.surv=TRUE)
  }
  
  The arguments are interpreted as
  \itemize{
    \item \code{pars} A vector of three parameters, in the order
    \code{b}, \code{mu}, \code{lambda}, which are speciation 
    initiation rate, extinction rate, speciation completion rate.
    \item \code{condition.surv} (logical): should the likelihood
    calculation condition on survival of two lineages and the speciation
    event subtending them?  This is done by default, following Nee et
    al. 1994.
  }
}

\examples{
pars <- c(0.3,0.1,0.1)
phy <- tree.prosse(pars,max.t=15)

##Sample one lineage per species
species.label <- unique(phy$species)
n.species <- length(species.label)
tip.species <- sapply(species.label,function(i) is.element(phy$species,i))
n.tip.species <- colSums(tip.species)
keep.tips <- sapply(1:n.species,function (i) ifelse(n.tip.species[i]==1,which(tip.species[,i]),sample(which(tip.species[,i]),1)))
to.drop <- rep(TRUE,length(phy$tip.label))
to.drop[keep.tips] <- FALSE
phy <- prune.prosse(phy,to.drop)

lik <- make.prosse.sp(phy)
p <- starting.point.prosse(phy)
fit <- find.mle(lik,p,lower=c(0,0,0))
}

\references{
  Hua X., Herdha T., and Burden C. 2021. Protracted speciation under the
  state-dependent speciation and extinction approach (bioRxiv).
}

\author{Xia Hua}
\keyword{models}

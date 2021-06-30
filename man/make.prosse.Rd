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
make.prosse(tree)
starting.point.prosse(tree)
}

\arguments{
  \item{tree}{An ultrametric bifurcating phylogenetic tree, in
    \code{ape} \dQuote{phylo} format. In the tree, tree$species, 
    named by tree$tip.label, gives the species identity of each tip}
}

\details{
  \code{make.prosse} returns a function of class \code{prosse}.

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

  \code{starting.point.bisse} produces a heuristic starting point to
  start from, based on the character-independent birth-death model. 
  Birth rate gives the starting point of speciation initiation rate
  and speciation completion rate. Death rate gives the starting point
  of extinction rate.
}

\examples{
pars <- c(0.3,0.1,0.1)
phy <- tree.prosse(pars,max.t=15)
lik <- make.prosse(phy)
p <- starting.point.prosse(phy)
fit <- find.mle(lik,p,lower=c(0,0,0))
}

\references{
  Hua X., Herdha T., and Burden C. 2021. Protracted speciation under the
  state-dependent speciation and extinction approach (bioRxiv).
}

\author{Xia Hua}
\keyword{models}

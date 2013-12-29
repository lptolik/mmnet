\name{stripSSN}
\alias{stripSSN}
\title{
strip the reference metabolic network
}
\description{
Given a set of KOs of a sample with a specific state and the reference metabolic network, this function creates and returns sub-network with only the supplied nodes (KOs) and any edges between them.
}
\usage{
stripSSN(refNet, abundance)
}
\arguments{
  \item{refNet}{
a igpah, the reference network saved in the \code{\link{RefDbcache}},
}
  \item{abundance}{
a numeric vector, abundance of KOs of a given specific state, can be obtained with \code{\link{estimateKOAbundance}}
}
}
\details{
Reference metabolic network has been saved in the RefDbcache. There are varies sets of KOs and its abundance between metagenomic samples. With this function, users can create a specified sub-network of the reference network for each sample.
}
\value{
A new graph object names mmnet, and the abundance is added as a vertex meta-attributes.
}

\seealso{
\code{\link{RefDbcache}}, \code{\link{estimateKOAbundance}}
}
\examples{
# load reference metabolic network 
data(RefDbcache)

# KOAbund
data(mgp10Anno)
KOAbund <- estimateKOAbundance(mgp10Anno[[1]])

subnet  <- stripSSN(RefDbcache$network,KOAbund)

}


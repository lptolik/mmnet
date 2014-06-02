\name{stripSSN}
\alias{stripSSN}
\title{
strip the reference metabolic network
}
\description{
Given a set of KOs of a sample with a specific state and the reference metabolic network, this function creates and returns sub-network with only the supplied nodes (KOs) and any edges between them.
}
\usage{
stripSSN(abundance)
}
\arguments{
\item{abundance}{ 
abundance of KOs of a given specific state, may a numeric vector which can be obtained with \code{\link{estimateAbundance}} 

or a given Biological Observation Matrix (BIOM) format, the ncol of BIOM file must be 1. The BIOM file format (canonically pronounced biome) is designed to be a general-use format for representing biological sample by observation contingency tables. More details on BIOM file can be found in \url{http://biom-format.org/}  R package \emph{biom}. There is also a example BIOM file in directory "inst/extdata" of \emph{mmnet}. 
}
}
\details{
Reference metabolic network has been saved in the RefDbcache. There are varies sets of KOs and its abundance between metagenomic samples. With this function, users can create a specified sub-network of the reference network for each sample.
}
\value{
A new graph object names mmnet will be return, and the abundance is added as a vertex meta-attributes.
}

\seealso{
\code{\link{RefDbcache}}, \code{\link{estimateAbundance}}
}
\examples{
# load reference metabolic network 
data(RefDbcache)

# take a numeric vector as input 
data(mgp10Anno)
abund <- estimateAbundance(mgp10Anno[[1]])
subnet  <- stripSSN(abund)

# BIOM format as input
if (require(biom)){
  abund2 <- read_biom(system.file("extdata", "enzymeAbund.biom", package = "mmnet"))
  subnet <- stripSSN(abund2)
}
}


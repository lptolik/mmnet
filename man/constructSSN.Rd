\name{constructSSN}
\alias{constructSSN}
\title{
Metabolic network construction of a specific sample
}
\description{
Once KO is estimated for a sample, a specific metabolic network can be constructed by calling this function.
}
\usage{
constructSSN(KOAbund, RefDbcache)
}
\arguments{
  \item{KOAbund}{
a vector, can be obationed by function \emph{estimateKOAbund}
}
\item{RefDbcache}{
reference data
}
}
\details{
A global network has been constructed first, more detalis in \code{\link{updateKEGGPathway}} and \code{link{constructMetabolicNetwork}}. \emph{constructSSN} constructs a specific network for each sample, users can compare arbitrary metabolic network. 
}
\value{
igraph, abundance as a node attributes
}

\seealso{
\code{\link{updateKEGGPathway}}, \code{link{constructMetabolicNetwork}}
}
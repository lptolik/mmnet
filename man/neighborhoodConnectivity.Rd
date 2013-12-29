\name{neighborhoodConnectivity}
\alias{neighborhoodConnectivity}
\title{
Neighborhood connectivity of vertices
}
\description{
The neighborhood connectivity of a node is defined as the average connectivity of all neighbors if it.  
}
\usage{
neighborhoodConnectivity(g, mode = c("in", "out", "all"))
}
\arguments{
  \item{g}{
the graph  to analyzed
}
  \item{mode}{
Character string, "out" for out-degree, "in" for in-degree or "all" for the sum of the two.
}
}

\value{
a numeric vector of the same length as number of vertices.
}


\seealso{
\code{\link{analyzeSSN}}
}
\examples{
##
  require(igraph)
  g <- graph.ring(10)
  neighborhoodConnectivity(g)
}

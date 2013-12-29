\name{neighborhoodConDst}
\alias{neighborhoodConDst}

\title{
Neighborhood connectivity distributin of vertices
}
\description{
Neighborhood connectivity distributin gives the average of the neighborhood connectivities of all vertex  with k neighbors.
}
\usage{
neighborhoodConDst(g, mode = c("in", "out", "all"))
}
\arguments{
   \item{g}{
the graph  to analyzed
}
  \item{mode}{
Character string, "out" for out-degree, "in" for in-degree or "all" for the sum of the two.
}
}


\seealso{
\code{\link{analyzeSSN}}
}
\examples{
  require(igraph)
  g <- graph.ring(10)
  neighborhoodConDst(g)
}


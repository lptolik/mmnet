\name{topologicalCoff}
\alias{topologicalCoff}
\title{
Topological coffecient of vertices 
}
\description{
Topological coffecient represents a relative measure for the extent to which a node shares neighbors with other nodes.
}
\usage{
topologicalCoff(g, mode = c("in", "out", "all"))
}
\arguments{
  \item{g}{
The graph to analyze.
}
  \item{mode}{
Character string, "out" for out-degree, "in" for in-degree or "all" for the sum of the two.
}
}
\details{
The topological coefficient of a vertex n is computed as followns: first calculte average neighbors of all vertex m that share at least one neighbor with n, then divided by the number of m. Moreover,  plus one if there is a direct link between n and m 
}
\value{
 a numeric vector of the same length as number of vertices.
}

\seealso{
\code{\link{analyzeSSN}}
}
\examples{
  require(igraph)
  g <- graph.ring(10)
  topologicalCoff(g)
}


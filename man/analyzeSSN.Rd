\name{analyzeSSN}
\alias{analyzeSSN}
\title{
Analyze the specific state network 
}
\description{
Calculate Topological features of each node and linked with abundance attributes.
}
\usage{
analyzeSSN(g, Scatterplot=TRUE, mode = c("in", "out", "all"), .properties = c(), ...)
}
\arguments{
  \item{g}{
a igraph network of a specific state, it can be obtained with function \code{\link{stripSSN}} 
}
  \item{Scatterplot}{
logical, if show the scatterplot between topolagical features and abundance
}
  \item{mode}{
Character constatnt, it specifies how to use the direction of the edges if a directed graph is analyzed. For "out" only the outgoing edges are followed, so all vertices reachable from the source vertex in at most order steps are counted. For "in" all vertices from which the source vertex is reachable in at most order steps are counted. "all" ignores the direction of the edges. This argument is ignored for undirected graphs.
}
  \item{.properties}{
nstead of specifying the topological features in free form via the ... argument, one can specify them as a vector
}
  \item{\dots}{
topological features that user want to calculate.
}
}
\details{
Five most presented vertex topological feartures is provided by this function.

betweenness centrality:  defined by the number of geodesics (shortest paths) going through a vertex , more details in\code{\link{betweenness}} of package igraph

closeness centrality:  measures how many steps is required to access every other vertex from a given vertex. more details in\code{\link{closeness}} of package igraph

degree: the most basic structural property, represents the number of its adjacent edges.

neighborhoodCon: Neighborhood Connectivity, is defined as the average connectivity of all neighbors of a vertex

topoparam:  Topological Coffecient is a relative measure for the extent to which a node 
shares neighbors with other nodes

}
\value{
a graph that add all the calculated topological features as the meta-attributes
}
\seealso{
\code{\link{stripSSN}}, \code{\link{topologicalCoff}},\code{link{neighborhoodConnectivity}}
}
\examples{
##
data(mgp10Anno)
data(RefDbcache)
KOAbund <- estimateKOAbundance(mgp10Anno[[1]])

subnet  <- stripSSN(RefDbcache$network,KOAbund)
analyzeSSN(subnet, Scatterplot=TRUE, .properties=list("betweenness","closeness","degree","    
          neighborhoodCon", "topoparam"))
}

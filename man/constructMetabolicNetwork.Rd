\name{constructMetabolicNetwork}
\alias{constructMetabolicNetwork}
\title{
Metabolic network construction of microbiome
}
\description{
Metabolic network construction according to the metabolites among enzymes
}
\usage{
constructMetabolicNetwork(RefDbcache)
}
\arguments{
  \item{RefDbcache}{
a list, contains the KO, its metabolites, reaction in all KEGG metabolic pathway as well as an initial reference network  
}
}

\value{
referce metabolic network of igraph
}


\examples{
## may take a long time
# data(RefDbcache)
# constructMetabolicNetwork(RefDbcache)
}

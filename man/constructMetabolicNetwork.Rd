\name{constructMetabolicNetwork}
\alias{constructMetabolicNetwork}
\title{
Construct a community-level metabolic network with the entire set of enzymes get from the sample
}
\description{
The function construct a community-level metabolic network between the enzymes. Nodes represent enzymes and edges connect enzymes with metabolites. 
}
\usage{
constructMetabolicNetwork(KEGGPathwayInfo, MGRAST.KEGGAnno,degree.cutoff=0)
}
\arguments{
  \item{KEGGPathwayInfo}{
The reference annotation database, includes information of enzymes and its metabolites from the KEGG metabolic pathway data. The KEGG database annotated enzymes with metabolic reactions, users will get the reference annotation database with function \code{getKEGGPathwayInfo}().
}
  \item{MGRAST.KEGGAnno}{
The enzyme (EC numbers) of the KEGG ontology annotation.
}
  \item{degree.cutoff}{an interger, nodes which degree is greater than it will be plot. The argument deletes the nodes and edges with no biological significance.
}

}
\details{
Users can construct a metabolic network with a set of enzymes according to the KEGG database. The KEGG database was used to annotate enzymes with metabolic reactions. Each enzymes may be associated with multiple enzymes. Using this mapping,
an enzyme-based metabolic network was constructed, where nodes represent enzymes and a directed edge from enzyme 1 to enzyme 2 indicates that a product metabolite of a reaction catalyzed by enzyme 1 is a substrate metabolite of a reaction catalyzed by enzyme 2.
}
\value{
Return a igraph graphs object.
}
\seealso{
\code{\link{getKEGGPathwayInfo}}, \code{\link{parseMGRASTAnno}}
}
\examples{
#load the reference database and example data
data(KEGGPathwayInfo)
data(test.ontology)

#ontology file parsing and network constructing
KEGGAnno=parseMGRASTAnno(test.ontology)
g=constructMetabolicNetwork(KEGGPathwayInfo,KEGGAnno,degree.cutoff=3)
plot.igraph(g,vertex.size=3,vertex.label=NA,
      	    edge.width=0.5,edge.arrow.size=0.1,
            vertex.color="green",
            edge.arrow.width=0.1,	
            layout=layout.fruchterman.reingold)

}


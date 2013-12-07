\name{KEGGPathwayInfo}
\alias{KEGGPathwayInfo}
\title{
Reference database (KEGG pathway) data information for constructing metabolic network
}
\description{
The data includes information of all enzymes, products and substrates from the KGML files of KEGG metabolic pathway with \code{getKEGGPathwayInfo}().
}
\usage{data(KEGGPathwayInfo)}
\format{
A list of length 3:
First, ec number, represents the ec number of enzymes in KEGG database.
Second, substrate, represents the corresponding substrates of enzymes.
Third, product, represents the corresponding products of enzymes.
}
\details{
The KEGG metabolic pathway is the best organized part of the KEGG/PATHWAY database, each metabolic pathway can be viewed as a network of enzymes or network of EC numbers. The KEGG Markup Language (KGML) is an exchange format of the KEGG pathway maps, which is converted from internally used KGML+ (KGML+SVG) format.It will retrieves all enzymes, products and substrates from the KGML files of KEGG metabolic pathway with \code{getKEGGPathwayInfo}(). This package saves the an initial enzymes and metabolites data in the subdirectory data named \code{\link{KEGGPathwayInfo}}. However, KEGG metabolic pathway updated frequently, users can update the file with function\code{\link{getKEGGPathwayInfo}}.
}
\references{
url{http://www.genome.jp/kegg/}
}
\examples{
data(KEGGPathwayInfo)
}
\keyword{datasets}

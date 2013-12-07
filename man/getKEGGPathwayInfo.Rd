\name{getKEGGPathwayInfo}
\alias{getKEGGPathwayInfo}
\title{
Updating the information of enzymes and its metabolites from the KEGG database 
}
\description{
Getting all the enzymes, the corresponding substrates and products involved in the KEGG metabolic 
pathway with the KEGG API form KEGG database.   
}
\usage{
getKEGGPathwayInfo()
}
\details{
The KEGG metabolic pathway is the best organized part of the KEGG/PATHWAY database, each metabolic pathway 
can be viewed as a network of enzymes or network of EC numbers. The KEGG Markup Language (KGML) is an exchange 
format of the KEGG pathway maps, which is converted from internally used KGML+ (KGML+SVG) format. It 
will retrieves all enzymes, products and substrates from the KGML files of KEGG metabolic pathway. This package 
saves the an initial enzymes and metabolites data in the subdirectory data named \code{\link{KEGGPathwayInfo}}.
 However, KEGG metabolic pathway updated frequently, the function will update the file 
\code{\link{KEGGPathwayInfo}}. It can be used as reference database of metabolic network constructing. 
For more details, see \code{\link{constructMetabolicNetwork}}.
}
\references{
\url{http://www.kegg.jp/kegg/xml/docs/}
}
\seealso{
\code{\link{constructMetabolicNetwork}},\code{\link{KEGGPathwayInfo}}
}
\examples{
getKEGGPathwayInfo()

}


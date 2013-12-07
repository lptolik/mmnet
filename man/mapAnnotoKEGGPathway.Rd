\name{mapAnnotoKEGGPathway}
\alias{mapAnnotoKEGGPathway}
\title{
Map a set of enzymes to the KEGG global metabolic pathway ec00010
}
\description{
The function map a set of interested enzymes to the KEGG global metabolic pathway \sQuote{\code{ec01100}} with the KEGG WebLinks. A new colored global metabolic pathway is will be opened in your broswer.
}
\usage{
mapAnnotoKEGGPathway(KEGGAnno, map.color = "blue",browse=FALSE)
}
\arguments{
  \item{KEGGAnno}{
A vector consists of EC numbers to be colored.
}
  \item{map.color}{
The map color of the enzyme in the KEGG global map. If it is numeric then the current palette is used, see \link{palette}. If it is a character vector then it may either contain named colors or RGB specified colors with three or four bytes. All strings starting with   \sQuote{\code{#}} are assumed to be RGB color specifications. The default is blue.
}
 \item{browse}{
logical, if TRUE, browse the url of the colored map, the default is FALSE.  
}
}
\details{
The most KEGG metabolic pathway maps are still most widely used including the KGML (KEGG XML) version. The globalmap ec0010 which is created as an SVG file by manually combining about all existing metabolic maps. This global metabolism mapallows the user to view and compare the entire metabolism, such as by mapping metagenomics data or microarray data. This function is coloring the interested enzymes on the global metabolic pathway.
}
\value{
Return a character string for the url of colored map.
}
\references{
url{http://www.kegg.jp/kegg/rest/weblink.html}
}
\seealso{
\code{\link{constructMetabolicNetwork}}
}
\examples{
data(test.ontology)
KEGGAnno=parseMGRASTAnno(test.ontology)
colored.map=mapAnnotoKEGGPathway(KEGGAnno,map.color="blue")
}

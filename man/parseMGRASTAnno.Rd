\name{parseMGRASTAnno}
\alias{parseMGRASTAnno}
\title{
Parsing the KEGG enzyme information from the MGRAST functional annotation files
}
\description{
The function parses MGRAST functional annotation files to retrieve the enzyme (EC numbers) of the KEGG 
ontology annotation. 
}
\usage{
parseMGRASTAnno(MGRAST.ontology)
}
\arguments{
  \item{MGRAST.ontology}{
A MGRAST functional annotation file suffixed with 'expand.ontology' , includes md5 value of the 
database sequence hit followed by sequence or cluster ID, similarity information, and annotation information.
 Users can download it in the MGRAST website.
}
}
\details{
With the annotation file, some values are encoded as index numbers in the annotation files of the argument.
Column fields are as follows:

    1. hit id (md5sum) index number
    2. query id (either fasta ID or cluster ID)
    3. percentage identity
    4. alignment length
    5. e-value
    6. function annotation index number
    7. functional category (ontology) index number
    8. DB source index number

This function is just to get the KEGG ontology from the ontology annotation files. Users can construct 
metabolic networks by mapping this annotation to the selected reference database (KEGG PATHWAY).
\code{\link{constructMetabolicNetwork}} for more details.
}
\value{
Return a named character vector of enzymes consists of EC numbers.
}
\references{
\url{http://metagenomics.anl.gov/}
}
\seealso{
\code{\link{constructMetabolicNetwork}}, \code{\link{mapAnnotoKEGGPathway}}
}
\examples{
# get the test ontology annotation file
data(test.ontology)

#parse the annotation file
KEGGAnno=parseMGRASTAnno(test.ontology)
}

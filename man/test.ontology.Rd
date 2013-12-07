\name{test.ontology}
\alias{test.ontology}
\title{
An example MGRAST ontology data
}
\description{
A MGRAST functional annotation file suffixed with 'expand.ontology' , includes md5 value of the database sequence hit followed by sequence or cluster ID, similarity information, and annotation information. Users can download it in the MGRAST website.
}
\usage{data(test.ontology)}
\format{
A tab-dilimited file format. Column fields are as follows:

    1. hit id (md5sum) index number
    2. query id (either fasta ID or cluster ID)
    3. percentage identity
    4. alignment length
    5. e-value
    6. function annotation index number
    7. functional category (ontology) index number
    8. DB source index number
}
\source{
\url{http://metagenomics.anl.gov/}
}
\examples{
data(test.ontology)
}
\keyword{datasets}

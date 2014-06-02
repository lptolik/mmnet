\name{getMgrastAnnotation}
\alias{getMgrastAnnotation}
\title{
Get the annotation profile from MGRAST
}
\description{
The function downloads the annotation profile in MGRAST
}
\usage{
getMgrastAnnotation(MetagenomeID, evalue = 5, identity = 60, length = 15, resource = c(source = "KO", type = "ontology"), webkey = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MetagenomeID}{
character, the metagenome ID in mgrast
}
  \item{evalue}{
expectation negative exponent value threshold for saving hits, default is 5
}
  \item{identity}{
minimum percent identity of matches, default is 60  
}
  \item{length}{
value for minimum alignment length cutoff, default is 15  
}
  \item{resource}{
a two length vecter, the source you want to download, e.g., 
c("KO","ontology") represents the ontology annotation with the KO database
, more details in MGRAST
}
  \item{webkey}{
private metagenomic data in MGRAST need a webkey to download
}
}
\value{
annotation profile
}
\references{
\url{http://api.metagenomics.anl.gov/1/api.html#annotation}
}

\examples{
##
# getMgrastAnnotation("mgm4447943.3")
}

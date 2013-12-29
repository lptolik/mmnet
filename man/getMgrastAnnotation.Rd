\name{getMgrastAnnotation}
\alias{getMgrastAnnotation}
\title{
Get the annotation profile from MGRAST
}
\description{
The function downloads the annotation profile in MGRAST
}
\usage{
getMgrastAnnotation(MetagenomeID, resource = c(source = "KO", type = "ontology"), webkey = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MetagenomeID}{
character, the metagenome ID in mgrast
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

\examples{
##
# getMgrastAnnotation("mgm4447943.3", resource=c("KO","ontology"), webkey = NULL)
}

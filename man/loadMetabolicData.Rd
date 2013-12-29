\name{loadMetabolicData}
\alias{loadMetabolicData}
\title{
Metabolic data loading
}
\description{
This function load the updated reference metabolic data into workspace.
}
\usage{
loadMetabolicData(path = Sys.getenv("HOME"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
 character, the path that reference data saved, default is the user's HOME directory.
}
}
\seealso{\code{link{updateKEGGPathway}}
}
\examples{
##
# loadMetabolicData(path = Sys.getenv("HOME"))
}
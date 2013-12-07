\name{mmnet-package}
\alias{mmnet-package}
\title{
Introduction to the mmnet Package
}
\description{
The \pkg{mmnet} add-on is an implementation of mapping metagenomic data to the KEGG global metabolic pathway and constructing a systems-level network. 
It is of great help to analysis the various functional properties, including regulation and metabolic functionality of the metagenome.
} 
\details{
This package contains functions for mapping metagenomic data to the KEGG global metabolic pathway and constructing systems-level network. First,  the information of enzymes, corresponding substrates and products is parsed from the KEGG PATYWAY database with function \code{\link{getKEGGPathwayInfo}}. Second, functional ontology annotation with KEGG is got by parse the MGRAST ontology annotation files with function \code{\link{parseMGRASTAnno}}. Finally, colored global metabolic pathway and systems-level network are made to analyze the metabolic functionality of metagenome with function \code{\link{mapAnnotoKEGGPathway}} and \code{\link{constructMetabolicNetwork}}.
}

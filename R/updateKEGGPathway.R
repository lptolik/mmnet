updateKEGGPathway <- function(path = Sys.getenv("HOME")) {
    ko.list <- keggList("pathway/ko")
    deleteGlobal.info <- regexpr("ko01|ko00", names(ko.list), perl = TRUE)
    ko.list <- ko.list[deleteGlobal.info > 0|(grepl('(B|b)acter',ko.list)&!grepl('(E|e)pithel',ko.list))]
    ko.path <- regexpr("ko\\d+", names(ko.list), perl = TRUE)
    ko.pathway <- substring(names(ko.list), ko.path, ko.path + attr(ko.path, "match.length") -
        1)
    pathway.info <- lapply(ko.pathway, getKOPathwayInfo)
    names(pathway.info) <- ko.pathway
    RefDbcache <- as.list(preprocessKOMetabolites(pathway.info))
    RefDbcache$pathway<-pathway.info
    ko.df<-data.frame(ko=gsub('^path:','',names(ko.list)),name=ko.list)
    RefDbcache$path.names<-ko.df
    saveMetabolicData(RefDbcache, path)
}

stripSSN <- function(abundance) {
    if (!exists("RefDbcache"))
        data(RefDbcache)
    RefDbcache <- get("RefDbcache", envir = parent.frame())
    if (!is.igraph(RefDbcache$network)) 
        stop("not a igraph object")
    if (class(abundance) != "numeric" || !is.vector(abundance)) 
        stop("abundance must be numeric of vector")
    refnode <- V(RefDbcache$network)$name
    subnodes <- intersect(names(abundance), refnode)
    if (!length(subnodes)) 
        stop("names of abundance should be KO number or there is no KO intersection between this sample and reference data")
    g <- induced.subgraph(RefDbcache$network, subnodes)
    match.index <- match(V(g)$name, names(abundance))
    subabund <- abundance[match.index]
    g <- set.graph.attribute(g, "name", "SSN")
    g <- set.vertex.attribute(g, "abundance", index = V(g), value = subabund)
    g <- delete.vertices(g, names(which(igraph::degree(g, mode = "all") == 0)))
    return(g)
} 

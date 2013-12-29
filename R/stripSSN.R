stripSSN <- function(refNet, abundance) {
    if (!is.igraph(refNet)) 
        stop("not a igraph object")
    if (class(abundance) != "numeric" || !is.vector(abundance)) 
        stop("abundance must be numeric of vector")
    refnode <- V(refNet)$name
    subnodes <- intersect(names(abundance), refnode)
    if (!length(subnodes)) 
        stop("names of abundance should be KO number or there is no KO intersection between this sample and reference data")
    g <- induced.subgraph(refNet, subnodes)
    match.index <- match(V(g)$name, names(abundance))
    subabund <- abundance[match.index]
    g <- set.graph.attribute(g, "name", "subNet")
    g <- set.vertex.attribute(g, "abundance", index = V(g), value = subabund)
    g <- delete.vertices(g, names(which(igraph::degree(g) == 0)))
    return(g)
} 

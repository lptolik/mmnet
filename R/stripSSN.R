stripSSN <- function(abundance) {
    if (!exists("RefDbcache"))
        data(RefDbcache)
    RefDbcache <- get("RefDbcache", envir = parent.frame())
    if (!is.igraph(RefDbcache$network)) 
        stop("not a igraph object")
    if (!inherits(abundance, "biom") && !is.vector(abundance)) 
        stop("abundance must be numeric of vector or BIOM format")
    refnode <- V(RefDbcache$network)$name
    if (class(abundance) == "biom"){
      if (ncol(abundance) != 1)
        stop("the ncol of biom file must be 1")
      if(abundance$matrix_type == "sparse"){
        abundance = as.matrix(biom_data(abundance))
      }else{
        abundance = biom_data(abundance) 
      }  
      subnodes <- intersect(rownames(abundance)[abundance != 0], refnode)
      g <- induced.subgraph(RefDbcache$network, subnodes)
      match.index <- match(V(g)$name, rownames(abundance))
    }else{
      subnodes <- intersect(names(abundance)[abundance != 0], refnode)
      g <- induced.subgraph(RefDbcache$network, subnodes)
      match.index <- match(V(g)$name, names(abundance))
    }
    if (!length(subnodes)) 
        stop("names of abundance should be KO number or there is no KO intersection between this sample and reference data")
    subabund <- abundance[match.index]
    g <- set.graph.attribute(g, "name", "SSN")
    g <- set.vertex.attribute(g, "abundance", index = V(g), value = subabund)
    g <- delete.vertices(g, names(which(igraph::degree(g, mode = "all") == 0)))
    return(g)
} 

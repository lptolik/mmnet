constructSSN <- function(KOAbund, RefDbcache) {
    ko <- RefDbcache$ko
    product <- RefDbcache$product
    substrate <- RefDbcache$substrate
    message("Constructing nodes of the network  ...", domain = NA)
    node1 <- intersect(ko, names(KOAbund))
    node1.index <- mapply(function(node1) which(ko == node1), node1)
    sub.substrate <- substrate[node1.index]
    sub.product <- product[node1.index]
    message("Constructing edges of the network  ...", domain = NA)
    edge.matrix <- matrix(0, length(sub.product), length(sub.product))
    edge.matrix <- mapply(function(y) mapply(function(x) intersect(y, x), sub.substrate), 
        sub.product)
    edge.matrix[edge.matrix > 1] <- 1
    rownames(edge.matrix) <- colnames(edge.matrix) <- node1
    diag(edge.matrix) <- 0
    message("constrcuting the ref Network ")
    g <- graph.adjacency(edge.matrix, mode = "directed")
    g <- set.graph.attribute(g, name = "name", value = "SSN")
    return(g)
} 

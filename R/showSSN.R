showSSN <- function(SSN, mode = c("ref", "strip", "compared"), method = c("OR", "rank", 
    "JSD"), cutoff, ...) {
    ## vertex color for compared network
    VertexColor <- function(value, method, cutoff) {
        vertex.color <- rep("grey", length(value))
        if (method == "OR") {
            value <- abs(log2(value))
            enrich.no <- which(value > cutoff[[2]])
            deplete.no <- which(value < cutoff[[1]])
        } else {
            deplete.no <- which(value < quantile(value, probs = cutoff[[1]]))
            enrich.no <- which(value > quantile(value, probs = cutoff[[2]]))
        }
        vertex.color[enrich.no] <- "red"
        vertex.color[deplete.no] <- "green"
        return(vertex.color)
    }
    
    SSN <- delete.vertices(SSN, names(which(igraph::degree(SSN) == 0)))
    method <- match.arg(method, c("OR", "rank", "JSD"))
    if (!is.igraph(SSN)) 
        stop("not a igraph object")
    mode <- match.arg(mode, c("ref", "strip", "compared"))
    if (mode == "strip") {
        tmp <- get.vertex.attribute(SSN, name = "abundance", index = V(SSN)) + 1
        if (is.null(tmp)) 
            stop("error mode, please check your graph mode")
        vertex.size2 <- c(1:length(tmp))
        vertex.size2[which(tmp == 1)] <- 2
        vertex.size2[which(tmp != 1)] <- 2 + log(tmp[which(tmp != 1)])
        plot.igraph(SSN, vertex.size = vertex.size2, ...)
    }
    if (mode == "ref") {
        plot.igraph(SSN, ...)
    }
    if (mode == "compared") {
        if (method == "OR") {
            OR <- get.vertex.attribute(SSN, name = "OR", index = V(SSN))
            if (is.null(OR)) 
                stop("error mode, please check your graph mode")
            vertex.color <- VertexColor(OR, method = "OR", cutoff = c(0.5, 2))
            plot.igraph(SSN, vertex.color = vertex.color, ...)
        }
        if (method == "rank") {
            diffabund <- get.vertex.attribute(SSN, name = "diffabund", index = V(SSN))
            if (is.null(diffabund)) 
                stop("error mode, please check your graph mode")
            vertex.color <- VertexColor(diffabund, method = "rank", cutoff = c(0.1, 
                0.9))
            plot.igraph(SSN, vertex.size = 3, vertex.color = vertex.color, ...)
        }
        if (method == "JSD") {
            JSD <- get.vertex.attribute(SSN, name = "JSD", index = V(SSN))
            if (is.null(JSD)) 
                stop("error mode, please check your graph mode")
            vertex.color <- VertexColor(JSD, method = "JSD", cutoff = c(0.1, 0.9))
            if (is.null(JSD)) 
                stop("error mode, please check your graph mode")
            plot.igraph(SSN, vertex.color = vertex.color, ...)
        }
    }
    
} 

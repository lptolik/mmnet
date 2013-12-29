analyzeSSN <- function(g, Scatterplot = TRUE, mode = c("in", "out", "all"), .properties = c(), 
    ...) {
    if (!is.igraph(g)) {
        stop("Not a igraph object")
    }
    if (missing(.properties)) {
        .properties <- c(...)
    } else {
        .properties <- c(.properties, c(...))
        if (length(.properties) == 0) 
            stop("No inputs passed to properties")
    }
    proper.options <- c("betweenness", "closeness", "degree", "neighborhoodCon", 
        "topoparam")
    property <- match(proper.options, unlist(.properties), nomatch = 0)
    if (max(property) == 0) 
        stop("input is invalid", domain = NULL)
    mode <- match.arg(mode, c("in", "out", "all"))
    abundance <- get.vertex.attribute(g, name = "abundance", index = V(g))
    ftopo.value <- function(x, type) {
        switch(type, bc = betweenness(x, v = V(x), directed = TRUE, weights = NULL, 
            nobigint = TRUE, normalized = TRUE), cn = closeness(x, vids = V(x), mode = mode, 
            weights = NULL, normalized = FALSE), dg = igraph::degree(x, v = V(x), 
            mode = mode, loops = TRUE, normalized = FALSE), nc = neighborhoodConnectivity(x, 
            mode = mode), tp = topologicalCoff(x, mode = mode))
    }
    topo.params <- data.frame(matrix(0, length(abundance), length(.properties) + 
        1))
    colnames(topo.params) <- c("abundance", unlist(.properties))
    topo.params$abundance <- abundance
    property <- which(property != 0)
    if (match(1, property, nomatch = 0)) {
        bc <- ftopo.value(g, "bc")
        topo.params$betweenness <- bc
        g <- set.vertex.attribute(g, name = "betweenness", index = V(g), value = bc)
    }
    if (match(2, property, nomatch = 0)) {
        cn <- ftopo.value(g, "cn")
        topo.params$closeness <- cn
        g <- set.vertex.attribute(g, name = "closeness", index = V(g), value = cn)
    }
    if (match(3, property, nomatch = 0)) {
        dg <- ftopo.value(g, "dg")
        topo.params$degree <- dg
        g <- set.vertex.attribute(g, name = "degree", index = V(g), value = dg)
    }
    if (match(4, property, nomatch = 0)) {
        nc <- ftopo.value(g, "nc")
        topo.params$neighborhoodCon <- nc
        g <- set.vertex.attribute(g, name = "neighborhoodCon", index = V(g), value = nc)
    }
    if (match(5, property, nomatch = 0)) {
        tp <- ftopo.value(g, "tp")
        topo.params$topoparam <- tp
        g <- set.vertex.attribute(g, name = "topoparam", index = V(g), value = tp)
    }
    if (Scatterplot) {
        dat <- melt(topo.params, "abundance")
        print(ggplot(dat, aes_string(x = "abundance", y = "value")) + geom_point() + stat_smooth(method = "lm", 
            colour = "red") + facet_wrap(~variable, ncol = 2, scales = "free") + 
            scale_x_log10())
    }
    return(network = g)
} 

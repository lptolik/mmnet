compareSSN <- function(..., method = c("OR", "rank", "JSD"), cutoff, p.value = TRUE, 
                           Visualization = TRUE) {
    ## odds ratio
    odds.ratio <- function(abund) {
        fodds <- function(x){
            if (is.null(dim(x)))
                return(sapply(x, function(y)y/(sum(x) - y)))
            else
                return(apply(x, 1, function(y) sum(y)/(sum(x) - sum(y))))
        }
        ratio <- lapply(abund,fodds)
        return(odds.ratio = Reduce("/", ratio))
    }

    ## Jensen Shannon divergence from KL divergence
    JSdiv <- function(abund.rank) {
        M <- (Reduce("+", abund.rank))/2
        KLD <- lapply(abund.rank, function(rank) mapply(function(x, y) KLdiv(cbind(x, 
            y))[1, 2], data.frame(t(rank)), data.frame(t(M))))
        KLD <- rowMeans(do.call(cbind, KLD))
        
    }
    if (!exists("RefDbcache"))
	    data(RefDbcache)
    RefDbcache <- get("RefDbcache", envir = parent.frame())
    ko.abund <- c(...)
    state <- unique(names(ko.abund))
    method <- match.arg(method, c("OR", "rank", "JSD"))
    if (length(state) == 0) 
        stop("input should have names attributes: the state of the metagenome")
    if (length(state) != 2) 
        stop("state of the metagenome must be 2", domain = NULL)
    index <- lapply(c(1, 2), function(x) which(match(names(ko.abund), state) == x))
    names(index) <- state
    ## pool all samples together
    kos <- lapply(ko.abund, names)
    nodes <- unique(unlist(kos))
    diff.ko <- lapply(kos, function(x) setdiff(nodes, x))
    diff.abund <- lapply(diff.ko, function(x) {
        abund <- rep(0, length(x))
        names(abund) <- x
        return(abund)
    })
    extend.abund <- data.frame(mapply(function(x, y) c(x, y), ko.abund, diff.abund))
    extend.kos <- data.frame(mapply(function(x, y) c(x, y), kos[-1], diff.ko[-1]))
    kos.index <- data.frame(apply(extend.kos, 2, function(x) match(nodes, x)))
    extend.abund[-1] <- mapply(function(x, y) x[y], extend.abund[-1], kos.index)
    extend.abund <- data.frame((sapply(extend.abund, function(x) x/sum(x))))
    abund <- lapply(index, function(x) extend.abund[, x])
    nodes2 <- intersect(nodes, V(RefDbcache$network)$name)
    g <- induced.subgraph(RefDbcache$network, intersect(nodes2, V(RefDbcache$network)$name))
    if (p.value) {
        if (length(ko.abund) < 10) 
            warning("p value for small Size of metagenome make no sense")
        p.value <- mapply(function(x, y) wilcox.test(x, y)$p.value, data.frame(t(abund[[1]])), 
            data.frame(t(abund[[2]])))
        names(p.value) <- nodes
        p.value <- p.value[match(V(g)$name, nodes)]
        
        g <- set.vertex.attribute(g, name = "p.value", value = p.value, index = V(g))
    }
    if (method == "OR") {
        OR <- odds.ratio(abund)
        OR <- OR[match(V(g)$name, nodes)]
        g <- set.vertex.attribute(g, name = "OR", value = OR, index = V(g))
        if (Visualization) 
            showSSN(g, mode = "compared", method = "OR", cutoff = c(0.5, 2), vertex.label = NA, 
                edge.width = 0.3, edge.arrow.size = 0.1, edge.arrow.width = 0.1, 
                layout = layout.fruchterman.reingold, vertex.size = 3)
    }
    if (method == "rank") {
        abund.rank <- lapply(-extend.abund, rank, na.last = TRUE)
        abund.rank <- lapply(index, function(x) do.call(cbind, abund.rank[x]))
        mean.rank <- lapply(abund.rank, apply, 1, mean)
        diff.rank <- Reduce("-", mean.rank)
        diff.rank <- diff.rank[match(V(g)$name, nodes)]
        g <- set.vertex.attribute(g, name = "diffabund", value = diff.rank, index = V(g))
        if (Visualization) 
            showSSN(g, mode = "compared", method = "rank", cutoff = c(0.1, 0.9), 
                vertex.label = NA, edge.width = 0.3, edge.arrow.size = 0.1, edge.arrow.width = 0.1, 
                layout = layout.fruchterman.reingold, vertex.size = 3)
    }
    if (method == "JSD") {
        if (length(ko.abund) < 10) 
            warning("Jensen Shannon divergence for small size of metagenome make no sense")
        abund.rank <- lapply(extend.abund, rank, na.last = TRUE)
        abund.rank <- lapply(abund.rank, function(x) x/sum(x))
        abund.rank <- lapply(index, function(x) do.call(cbind, abund.rank[x]))
        JSD <- JSdiv(abund.rank)
        JSD <- JSD[match(V(g)$name, nodes)]
        g <- set.vertex.attribute(g, name = "JSD", value = JSD, index = V(g))
        if (Visualization) 
            showSSN(g, mode = "compared", method = "JSD", cutoff = c(0.1, 0.9), vertex.label = NA, 
                edge.width = 0.3, edge.arrow.size = 0.1, edge.arrow.width = 0.1, 
                layout = layout.fruchterman.reingold, vertex.size = 3)
    }
    return(g)
} 
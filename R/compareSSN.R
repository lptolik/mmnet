compareSSN <- function(abundance, method = c("OR", "rank", "JSD"), cutoff, p.value = TRUE, 
                           Visualization = TRUE) {
    ## odds ratio
    odds_ratio <- function(abund) {
        fodds <- function(x){
            if (is.null(dim(x)))
                return(sapply(x, function(y)y/(sum(x) - y)))
            else
                return(apply(x, 1, function(y) sum(y)/(sum(x) - sum(y))))
        }
        ratio <- lapply(abund,fodds)
        return(odds.ratio = Reduce("/", ratio))
    }
    if (!exists("RefDbcache"))
	    data(RefDbcache)
    RefDbcache <- get("RefDbcache", envir = parent.frame())
    ##ko.abund <- c(...)
    ko.abund <- data.frame(as.matrix(biom_data(abundance)))
    kos <- rownames(ko.abund)
    ko.abund <- data.frame(lapply(ko.abund,function(x)x/sum(x)))
    ## odds ratio
    
    colnames(ko.abund) <- gsub("\\.\\d+","", colnames(ko.abund))
    state <- unique(colnames(ko.abund))
    method <- match.arg(method, c("OR", "rank", "JSD"))
    if (length(state) == 0) 
        stop("input should have names attributes: the state of the metagenome")
    if (length(state) != 2) 
        stop("state of the metagenome must be 2", domain = NULL)
    index <- lapply(c(1, 2), function(x) which(match(colnames(ko.abund), state) == x))
    names(index) <- state
    g <- induced.subgraph(RefDbcache$network, intersect(kos, V(RefDbcache$network)$name))
    abund <- lapply(index, function(x) ko.abund[, x])
    if (p.value) {
        if (ncol(ko.abund) < 10) 
            warning("p value for small Size of metagenome make no sense")
        p.value <- mapply(function(x, y) wilcox.test(x, y)$p.value, data.frame(t(abund[[1]])), 
            data.frame(t(abund[[2]])))
        names(p.value) <- kos
        p.value <- p.value[match(V(g)$name, kos)]
        
        g <- set.vertex.attribute(g, name = "p.value", value = p.value, index = V(g))
    }
    if (method == "OR") {
        OR <- odds_ratio(abund)
        OR <- OR[match(V(g)$name, kos)]
        g <- set.vertex.attribute(g, name = "OR", value = OR, index = V(g))
        if (Visualization) 
            showSSN(g, mode = "compared", method = "OR", cutoff = c(0.5, 2), vertex.label = NA, 
                edge.width = 0.3, edge.arrow.size = 0.1, edge.arrow.width = 0.1, 
                layout = layout.fruchterman.reingold, vertex.size = 3)
    }
    if (method == "rank") {
        abund.rank <- lapply(ko.abund, order, na.last = TRUE, decreasing = TRUE)
        abund.rank <- lapply(abund.rank, order, na.last = TRUE)
        abund.rank <- lapply(index, function(x) do.call(cbind, abund.rank[x]))
        mean.rank <- lapply(abund.rank, apply, 1, mean)
        diff.rank <- Reduce("-", mean.rank)
        diff.rank <- diff.rank[match(V(g)$name, kos)]
        g <- set.vertex.attribute(g, name = "diffabund", value = diff.rank, index = V(g))
        if (Visualization) 
            showSSN(g, mode = "compared", method = "rank", cutoff = c(0.1, 0.9), 
                vertex.label = NA, edge.width = 0.3, edge.arrow.size = 0.1, edge.arrow.width = 0.1, 
                layout = layout.fruchterman.reingold, vertex.size = 3)
    }
    if (method == "JSD") {
        if (ncol(ko.abund) < 10) 
            warning("Jensen Shannon divergence for small size of metagenome make no sense")
        abund.rank <- lapply(ko.abund, order, na.last = TRUE, decreasing = TRUE)
        abund.rank <- lapply(abund.rank, order, na.last = TRUE)
        abund.rank <- lapply(abund.rank, function(x) x/sum(x))
        abund.rank <- lapply(index, function(x) do.call(cbind, abund.rank[x]))
        #JSD <- JSdiv(abund.rank)
        JSD <- mapply(function(x, y) sum(flexmix::KLdiv(cbind(x, y)))/2, data.frame(t(abund.rank[[1]]), row.names = NULL), 
                      data.frame(t(abund.rank[[2]]), row.names = NULL))
        JSD <- JSD[match(V(g)$name, kos)]
        g <- set.vertex.attribute(g, name = "JSD", value = JSD, index = V(g))
        if (Visualization) 
            showSSN(g, mode = "compared", method = "JSD", cutoff = c(0.1, 0.9), vertex.label = NA, 
                edge.width = 0.3, edge.arrow.size = 0.1, edge.arrow.width = 0.1, 
                layout = layout.fruchterman.reingold, vertex.size = 3)
    }
    return(g)
} 
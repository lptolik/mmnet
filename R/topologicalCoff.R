topologicalCoff <- function(g, mode = c("in", "out", "all")) {
    if (!is.igraph(g)) {
        stop("Not a igraph object")
    }
    mode <- match.arg(mode, c("in", "out", "all"))
    neigh <- neighborhood(g, 1, mode = mode)
    f <- function(z) {
        if (length(z) <= 2) {
            return(0)
        } else {
            node.neighbor <- z[-1]
            neigh1 <- sapply(node.neighbor, function(node.neighbor) neighborhood(g, 
                1, nodes = node.neighbor, mode = mode))
            del.self <- unlist(lapply(neigh1, function(x) x[-1]))
            subnode.neighbor <- as.data.frame(table(del.self), stringsAsFactors = FALSE)
            subnode.neighbor[, 1] <- as.numeric(subnode.neighbor[, 1])
            fneighborbind <- function(x) {
                index <- match(x, subnode.neighbor[, 1], nomatch = 0)
                if (index) {
                  subnode.neighbor[index, 2] <- subnode.neighbor[index, 2] + 1
                } else {
                  subnode.neighbor <- rbind(subnode.neighbor, c(x, 1))
                }
                return(subnode.neighbor)
            }
            for (i in 1:length(node.neighbor)) subnode.neighbor <- fneighborbind(node.neighbor[i])
            if (match(z[1], subnode.neighbor[, 1], nomatch = 0)) 
                subnode.neighbor <- subnode.neighbor[-match(z[1], subnode.neighbor[, 
                  1], nomatch = 0), ]
            return(sum(subnode.neighbor[, 2])/length(node.neighbor))
        }
    }
    topo.para <- sapply(neigh, f)
    return(topo.para)
} 

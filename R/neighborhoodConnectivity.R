neighborhoodConnectivity <- function(g, mode = c("in", "out", "all")) {
    if (class(g) != "igraph") 
        stop("g must be igraph")
    mode <- match.arg(mode, c("in", "out", "all"))
    neigh <- neighborhood(g, 1, mode = mode)
    fnc <- function(x) {
        if (length(x) >= 2) {
            neighSize <- lapply(x[-1], function(y) (neighborhood.size(g, 1, nodes = y, 
                mode = mode) - 1))
            return(mean(unlist(neighSize)))
        } else {
            return(0)
        }
    }
    return(sapply(neigh, fnc))
} 

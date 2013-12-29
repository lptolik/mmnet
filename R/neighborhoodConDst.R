neighborhoodConDst <- function(g, mode = c("in", "out", "all")) {
    if (!is.igraph(g)) {
        stop("Not a igraph object")
    }
    mode <- match.arg(mode, c("in", "out", "all"))
    neighCoff <- neighborhoodConnectivity(g, mode = mode)
    nbhsize <- neighborhood.size(g, 1, mode = mode) - 1
    df <- data.frame(nbhNo = nbhsize, nbhCoff = neighCoff)
    nbhCntDtb <- aggregate(df[, 2], list(df[, 1]), mean)
    names(nbhCntDtb) <- c("number.neighbor", "avg.neighborhoodCon")
    return(nbhCntDtb)
} 

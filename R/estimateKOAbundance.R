estimateKOAbundance <- function(KOAnno) {
    if (ncol(KOAnno) != 13) 
        stop("annotation file is invalid")
    if (grepl("query", KOAnno[1, 1])) {
        colnames(KOAnno) <- data.frame(lapply(KOAnno[1, ], as.character), stringsAsFactors = FALSE)
        KOAnno <- tail(KOAnno, -1)
    } else {
        stop("the first row should be the description of data")
    }
    if (grep("Download\\s+complete", KOAnno[nrow(KOAnno), 1])) {
        KOAnno <- head(KOAnno, -1)
    } else {
        stop("the last row should be the tag of data")
    }
    seq.ko <- KOAnno[, c(1, 13)]
    if (length(grep("^K", seq.ko[, 2])) != nrow(seq.ko)) 
        stop("all reads should be annotated with KO")
    ko.aggreg <- aggregate(seq.ko[, 2], list(seq.ko[, 1]), paste)
    ko <- sapply(sapply(sapply(ko.aggreg[, 2], function(x) strsplit(x, ";")), unlist), 
        unique)
    ko.aggreg[, 3] <- listLen(ko)
    ko.score <- data.frame(unlist(ko), rep(1/listLen(ko), listLen(ko)), stringsAsFactors = F)
    ko.count <- aggregate(ko.score[, 2], list(ko.score[, 1]), sum)
    ko.score[, 2] <- 1
    ko.reads <- aggregate(ko.score[, 2], list(ko.score[, 1]), sum)
    ko.abundance <- data.frame(cbind(ko.count, ko.reads[, 2]))
    names(ko.abundance) <- c("ko", "abundance", "reads")
    reads <- NULL
    df <- subset(ko.abundance, reads > 2)
    ret <- df[, 2]
    names(ret) <- df[, 1]
    return(ret)
} 

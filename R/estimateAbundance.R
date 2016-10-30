estimateAbundance <- function(KOAnno) {
    estimateSingleAnno <- function(Anno){
      if (ncol(Anno) != 13)
          stop("annotation file from MG-RAST is invalid")
      # that part is already checked in getMgrast
      # if (grepl("query", Anno[1, 1])) {
      #     colnames(Anno) <- data.frame(lapply(Anno[1, ], as.character), stringsAsFactors = FALSE)
      #     Anno <- tail(Anno, -1)
      # } else {
      #     stop("the first row should be the description of data")
      # }
      # that part is already checked in getMgrast
      # if (grep("Download\\s+complete", Anno[nrow(Anno), 1])) {
      #     Anno <- head(Anno, -1)
      # } else {
      #     stop("the last row should be the tag of data")
      # }
      seq.ko <- Anno[, c(1, 13)]
      if (length(which(grepl("^K", seq.ko[, 2]))) != nrow(seq.ko)){
        seq.ko[,2]<-sapply(seq.ko[,2],function(.x)
          paste(
            gsub('accession=\\[K([0-9]+)\\].*',
                 'K\\1',
                 unlist(str_split(.x,';'))),
            collapse = ';'))
      }#stop("all reads should be annotated with KO")
      #ko.aggreg <- aggregate(seq.ko[, 2], list(seq.ko[, 1]), paste)
      seq.ko.dt<-setDT(seq.ko)
      ko.aggreg.dt <- seq.ko.dt[,.(id = list(V13)),by=.(V1)]
      l.dt<-ko.aggreg.dt[, .(id)]
      ko.dt <- sapply(sapply(lapply(l.dt, function(x) str_split(x, ";")), unlist),
          unique)
      ko.aggreg.dt$len <- listLen(ko.dt)
      ko.score.dt <- data.frame(unlist(ko.dt), rep(1/listLen(ko.dt), listLen(ko.dt)), stringsAsFactors = F)
      ko.count.dt <- aggregate(ko.score.dt[, 2], list(ko.score.dt[, 1]), sum)
      ko.score.dt[, 2] <- 1
      ko.reads.dt <- aggregate(ko.score.dt[, 2], list(ko.score.dt[, 1]), sum)
      ko.abundance.dt <- data.frame(cbind(ko.count.dt, ko.reads.dt[, 2]))
      names(ko.abundance.dt) <- c("ko", "abundance", "reads")
      reads <- NULL
      df <- subset(ko.abundance.dt, reads > 2)
      ret <- df[, 2]
      names(ret) <- df[, 1]
      return(ret)
    }
    if (!is.data.frame(KOAnno)){
      ko.abund <- lapply(KOAnno, estimateSingleAnno)
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
      #extend.abund <- data.frame((sapply(extend.abund, function(x) x/sum(x))))
      #names(extend.abund) <- names(ko.abund)
      rownames(extend.abund) <- nodes
      #abund <- lapply(index, function(x) extend.abund[, x])
    }else{
      extend.abund <- estimateSingleAnno(KOAnno)
    }
    extend.abund <- data.frame(extend.abund)
    biom.data <- make_biom(extend.abund, observation_metadata = rownames(extend.abund))
    biom.data$type <- "enzymatic genes abundance"
    return(biom.data)
}

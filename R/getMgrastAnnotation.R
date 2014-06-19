getMgrastAnnotation <- function(MetagenomeID, evalue = 5, identity = 60, length = 15,
                                resource = c(source = "KO", type = "ontology"), login.info = NULL) 
{
    webkey <- login.info["webkey"]
    server.resource <- "http://api.metagenomics.anl.gov/1/annotation/similarity"
    server.para <- paste(paste0(c("?", rep("&", 4)), paste(c("source", "type", "evalue","identity", 
                                                            "length"), c(resource,evalue,identity,length), 
                                                          sep = "=")), collapse = "")
    url.str <- paste0(server.resource, "/", MetagenomeID, server.para)
    message("Loading the annotations form MG-RAST...", domain = NA)
    message("The time spent in this step is proportional to the total amount of remote data...")
    tryCatch(anno <- getURL(url.str),
        error = function(err) {stop(simpleError("Your Internet not connected or MGRAST host can not be connetecd, please try later"))}
    )
    private.index <- which(grepl("insufficient\\s+permissions", anno))
    if (length(private.index)) {
        if (is.null(login.info)) {
            stop("Load private data needs account login")
        } else {
            private.url <- paste0(url.str[private.index], "&auth=", webkey)
            anno[private.index] <- getURL(private.url)
        }
    }
    invalid.index <- which(grepl("does\\s+not\\s+exists", anno))
    invalid.source <- which(grepl("Invalid\\s+ontology\\s+source", anno))
    if (length(invalid.source)) 
        stop("invalid ontology source")
    if (length(invalid.index)) 
        stop(paste("The", paste(invalid.index, collapse = ","), "metagenomeID does not exists"))
    if (length(which(grepl("insufficient\\s+permissions", anno)))) 
        stop("invalid webkey")
    # anno <- lapply(as.list(anno), function(x) read.delim(textConnection(x), header = FALSE, 
    #     sep = "\t", stringsAsFactor = F))
    anno <- read.delim(textConnection(anno), header = FALSE, sep = "\t", stringsAsFactor = F)
    return(anno)
} 

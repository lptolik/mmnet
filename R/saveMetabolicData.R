saveMetabolicData <- function(RefDbcache, path = Sys.getenv("HOME")) {
    RefDbcache <- as.list(RefDbcache)
    RefDbcache$user <- Sys.getenv("USERNAME")
    RefDbcache$date <- date()
    RefDbcache$version <- R.Version()
    RefDbcache$date <- date()
    RefDbcache$network <- constructMetabolicNetwork(RefDbcache)
    message("saving data to the Specified dir...", domain = NA)
    path <- sprintf("%s/%s", path, ".mmnet")
    if (!exists(path)) 
        dir.create(path, showWarnings = F)
    save(RefDbcache, file = sprintf("%s/%s", path, "RefDbcache.rda"))
} 

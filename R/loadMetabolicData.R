loadMetabolicData <- function(path = Sys.getenv("HOME")) {
    path <- sprintf("%s/%s", path, ".mmnet")
    if (!exists(path)) {
        message(sprintf("No data found in %s, attempt to load default data.\nRecommand run updateKEGGPathway() to get the newest version of data.", 
            path), domain = NA)
        data(RefDbcache)
    } else load(file = sprintf("%s/%s", path, "RefDbcache.rda"))
} 

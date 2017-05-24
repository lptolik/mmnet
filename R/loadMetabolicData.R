loadMetabolicData <- function(path = Sys.getenv("HOME")) {
    path <- sprintf("%s/%s", path, ".mmnet")
    if (!dir.exists(path)) {
        message(" Attempt to load default data.\nRecommend run updateKEGGPathway() to get the latest version of data."
                , domain = NA)
        data(RefDbcache)
    } else{ 
      load(file = sprintf("%s/%s", path, "RefDbcache.rda"))
      assign('RefDbcache',RefDbcache,envir = globalenv())
    }
}

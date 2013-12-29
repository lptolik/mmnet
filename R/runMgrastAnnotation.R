runMgrastAnnotation <- function(fn, user, userpwd, new_project) {
    ## check project status, FALSE represents in processing, TRUE is completely
    projectStatus <- function(metagenome.id) {
        url.response <- getURL(paste0("http://api.metagenomics.anl.gov/1/metagenome/", 
            paste0("mgm", metagenome.id), "?verbosity=metadata", "&auth=", login.info$webkey))
        if (grepl("ERROR", url.response)) 
            return(FALSE) else return(TRUE)
    }
    
    if (!file.exists(fn)) {
        message("Can not find fn, please check if fn exsits!")
        return(FALSE)
    }
    login.info <- loginMgrast(user, userpwd)
    file <- listMgrastInbox(login.info)
    if (fn %in% file$files) {
        if (fn %in% names(file$locks)) {
            stop("File is in computing suequence status!")
            return(FALSE)
        }
        stop("File exists on MGRAST, remove it first!")
        delMgrastInbox(login.info, fn)
    }
    uploadMgrast(login.info, fn)
    file <- listMgrastInbox(login.info)
    while (fn %in% file$locks) {
        # stop('The fn is under statistics calculation')
        Sys.sleep(10)
        file <- listMgrastInbox(login.info)
    }
    if (missing(new_project)) 
        new_project <- as.character(Sys.time())
    metagenome <- submitMgrastJob(login.info, seqfile = basename(fn), new_project = new_project)
    metagenome.id <- str_extract_all(metagenome[[1]], "\\d+\\.\\d+")
    status <- projectStatus(metagenome.id)
    ## create txtprogressbar
    pb <- txtProgressBar(width = 40)
    textBar <- function(x = sort(runif(5))) {
        for (i in c(0, x, 1)) {
            Sys.sleep(0.5)
            setTxtProgressBar(pb, i)
        }
    }
    
    init.value <- seq(0, 1, by = 0.2)
    message("In annotation, may take time form days to weeks...")
    while (!status) {
        textBar(init.value)
        status <- projectStatus(metagenome.id)
    }
    close(pb)
    message("Annotation is complete...")
    message("Download the annotation profile from MGRAST server...")
    anno <- getMgrastAnnotation(paste0("mgm", metagenome.id), resource = c(source = "KO", 
        type = "ontology"), webkey = login.info$webkey)
    message("Done...")
    return(list(mgrastID = paste0("mgm", metagenome.id), KOAnnotation = anno))
} 

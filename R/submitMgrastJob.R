submitMgrastJob <- function(login.info, seqfile, new_project) {
    dereplication <- "dereplication"
    screening <- "h_sapiens_asm"
    filter_ln <- "filter_ln"
    deviation <- 2
    filter_ambig <- "filter_amibg"
    max_ambig <- 5
    para <- c(create_job = 1, table_perpage_0 = 10, seqfiles = seqfile, dereplication = dereplication, 
        screening = screening, filter_ln = filter_ln, deviation = deviation, filter_ambig = filter_ambig, 
        max_ambig = max_ambig, priorityOption = "never")
    if (FALSE) {
        if (missing(project)) {
        }
        if (missing(new_project)) {
            para <- c(para, new_project = as.character(Sys.time()))
        } else {
            para <- c(para, project = project)
        }
    }
    para <- c(para, new_project = new_project)
    response <- getForm("http://metagenomics.anl.gov/metagenomics.cgi?page=Upload", 
        .params = para, curl = login.info$curlhandle)
    response <- htmlParse(response)
    user <- getNodeSet(response, "//div[@id='user']")
    account <- sapply(user, xmlValue)
    account <- gsub("\\n+", "", account)
    metagenome.id <- getNodeSet(response, "//div[@class='modal-body']")
    metagenome <- sapply(metagenome.id, function(x) sapply(x["p"], xmlValue))
    metagenome <- unlist(try(metagenome[which(grepl("success", metagenome))], silent = T))
    names(metagenome) <- NULL
    if (length(metagenome)) {
        print(metagenome[1])
    } else {
        stop("job submmission failed")
    }
    return(list(MgrastID = metagenome.id, account = account))
} 

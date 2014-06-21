submitMgrastJob <- function(login.info, seqfile, new_project) {
    file <- listMgrastInbox(login.info)
    if (!match(seqfile, file$files, nomatch = 0))
        stop(paste0(seqfile, " is not in your Mgrast Inbox"))
    ## unpack the sequence file
    websession <- login.info$session
    if (grepl("(\\.gz$)|(\\.zip$)|(\\.tgz$)", seqfile))
    {
        message("Unpack the compressed File...",domain = NULL)
        unpackfile <- getForm("http://metagenomics.anl.gov/upload.cgi/user_inbox/?callback=1", 
            auth = websession, websession = "faction", faction = "unpack", unpack = "fn", fn=seqfile)
        seqfile <- file_path_sans_ext(seqfile)
    }

    check.dup <- getForm("http://metagenomics.anl.gov/metagenomics.cgi?page=Upload",auth=login.info$session,
        action="check_for_duplicates",curl = login.info$curlhandle,
        seqfiles=seqfile)
    if (grepl("already exist in MG-RAST",check.dup)){
        warning(paste(seqfile, "already exist in MG-RAST","\n", "existing ID:", str_extract(check.dup, "\\d+\\.\\d+")))
        return(str_extract(check.dup, "\\d+\\.\\d+"))
    }
    file <- listMgrastInbox(login.info)
    message("Under statistics caculation...", domain = NULL)
    while (seqfile %in% names(file$locks)) {
    # stop('The fn is under statistics calculation')
    refresh <- getForm("http://metagenomics.anl.gov/upload.cgi/user_inbox/?callback=1", 
            auth = websession, websession = "faction")
    Sys.sleep(10)
    file <- listMgrastInbox(login.info)
   }
    message("Statistics caculation is complete...", domain = NULL)
    dereplication <- "dereplication"
    screening <- "h_sapiens_asm"
    filter_ln <- "filter_ln"
    deviation <- 2
    filter_ambig <- "filter_amibg"
    max_ambig <- 5
    # seqfile <- file_path_sans_ext(seqfile)
    if (missing(new_project)) 
        new_project <- as.character(Sys.time())
    para <- c(create_job = 1, table_perpage_0 = 10, seqfiles = seqfile, dereplication = dereplication, 
        screening = screening, filter_ln = filter_ln, deviation = deviation, filter_ambig = filter_ambig, 
        max_ambig = max_ambig, priorityOption = "never")
    para <- c(para, new_project = new_project)
    message("Create the MGRAST project...", domain= NA )
    tryCatch(response <- getForm("http://metagenomics.anl.gov/metagenomics.cgi?page=Upload", 
        .params = para, curl = login.info$curlhandle),
        error = function(err) {stop(simpleError("Your Internet not connected or MGRAST host can not be connetecd, please try later"))}
    )
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
    return(list(MgrastID = metagenome[2], account = account))
    message("Job submiited", domain= NA )
} 

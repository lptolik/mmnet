listMgrastInbox <- function(login.info) {
    if (FALSE) {
        webkey <- login.info$webkey
        file.inform <- getForm("http://api.metagenomics.anl.gov/1/inbox", auth = webkey)
        if (length(grep("unauthorized\\s+request", file.inform)) != 0) {
            stop("the authorize may incorrect")
        }
        return(fromJSON(file.inform))
    }
    websession <- login.info$session
    file.inform <- getForm("http://metagenomics.anl.gov/upload.cgi/user_inbox/?callback=1", 
        auth = websession)
    if (length(grep("unauthorized\\s+request", file.inform)) != 0) {
        stop("the authorize may incorrect")
    }
    str_sub(file.inform, -3, -1) <- ""
    str_sub(file.inform, end = 27) <- ""
    return(fromJSON(file.inform))
} 

listMgrastInbox <- function(login.info) {
    websession <- login.info$session
    tryCatch(
            file.inform <- getForm("http://metagenomics.anl.gov/upload.cgi/user_inbox/?callback=1", 
                    auth = websession),
             error = function(err) {stop(simpleError("Your Internet not connected or MGRAST host can not be connetecd, please try later"))}
    )
    if (length(grep("unauthorized\\s+request", file.inform)) != 0) {
        stop("the authorize may incorrect")
    }
    str_sub(file.inform, -3, -1) <- ""
    str_sub(file.inform, end = 27) <- ""
    return(fromJSON(file.inform))
} 

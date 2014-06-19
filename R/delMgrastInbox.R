delMgrastInbox <- function(login.info, sequence) {
    websession <- login.info$session
    tryCatch(delete.file <- getForm("http://metagenomics.anl.gov/upload.cgi/user_inbox/?callback=1", 
              auth = websession, websession = "faction", faction = "del", del = "fn", fn = sequence),
            error = function(err) {stop(simpleError("Your Internet not connected or MGRAST host can not be connetecd, please try later"))}
    )
    ret <- delete.file
    str_sub(ret, -3, -1) <- ""
    str_sub(ret, end = 27) <- ""
    ret <- fromJSON(ret)
    if (length(ret$popup_messages) > 0) 
         print(ret$popup_messages)
    return(!(sequence %in% names(ret$locks)))
} 

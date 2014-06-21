uploadMgrast <- function(login.info, file) {
    if (!file.exists(file)) {
        message("Can not find file, please check if file exsits!")
        return(FALSE)
    }
    if (any(grepl(file_path_sans_ext(basename(file)),listMgrastInbox(login.info)$files)))
        stop("sequence file have been exits in your Inbox")
    tryCatch( 
        web.feedback <- postForm("http://api.metagenomics.anl.gov/1/inbox/upload", 
            .params = list(upload = fileUpload(file)), .opts = list(httpheader = 
            c(auth = login.info$webkey),noprogress = FALSE)),
        error = function(err) {stop(simpleError("Your Internet not connected or MGRAST host can not be connetecd, please try later"))}
    )
    if (grepl("data\\s+received\\s+successfully", web.feedback)) {
        cat(paste0("\n", file, " successfully upload ", date(), "\n"))
        return(TRUE)
    } else {
        cat("please check your webkey and the file, and ensure uploading file does not exist in your inbox\n")
        return(FALSE)
    }
} 
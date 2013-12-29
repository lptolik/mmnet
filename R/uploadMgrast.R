uploadMgrast <- function(login.info, file) {
    web.feedback <- try(postForm("http://api.metagenomics.anl.gov/1/inbox/upload", 
        .params = list(upload = fileUpload(file)), .opts = list(httpheader = c(auth = login.info$webkey))), 
        silent = F)
    if (grepl("data\\s+received\\s+successfully", web.feedback)) {
        cat(paste0(file, " successfully upload ", date(), "\n"))
        return(TRUE)
    } else {
        cat("please check your webkey and the file, and ensure uploading file does not exist in your inbox\n")
        return(FALSE)
    }
} 

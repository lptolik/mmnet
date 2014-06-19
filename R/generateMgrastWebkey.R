generateMgrastWebkey <- function(login.info) {
    ch <- login.info
    if (class(ch) != "CURLHandle") 
        stop("login curlhandle is need for generate your webkey")
    tryCatch(
        webkey.response <- getForm("http://metagenomics.anl.gov/metagenomics.cgi?page=Upload", 
            action = "generate_webkey", generate_new_key = "1", curl = ch),
        error = function(err) {stop(simpleError("Your Internet not connected or MGRAST host can not be connetecd, please try later"))}
    )
    webkey.date <- regexpr("(\\d{4})\\s+(\\d{2})\\-(\\d{2})\\s+(\\d{2}):(\\d{2})\\.(\\d{2})", 
        webkey.response, , perl = T)
    webkey.date <- substring(webkey.response, webkey.date, webkey.date + attr(webkey.date, 
        "match.length") - 1)
    webkey <- regexpr("\\w{25}", webkey.response, perl = T)
    webkey <- substring(webkey.response, webkey, webkey + attr(webkey, "match.length") - 
        1)
    return(list(webkey = webkey, invalid.until = webkey.date))
} 

#' Title
#'
#' @param MetagenomeID
#' @param evalue
#' @param identity
#' @param length
#' @param resource
#' @param login.info
#' @param webkey
#'
#' @return
#' @export
#'
#' @examples
getMgrastAnnotation <- function(MetagenomeID, evalue = 5, identity = 60, length = 15,
                                resource = c(source = "KO", type = "ontology"), login.info = NULL,webkey=NULL)
{
  checkAnno <- function(Anno){
    if (length(str_split(ncol(Anno),'\t')) != 13)
      stop("annotation file from MG-RAST is invalid")
    if (grepl("query", Anno[1])) {
      cn<-data.frame(lapply(str_split(Anno[1],'\t'), as.character), stringsAsFactors = FALSE)
      Anno <- tail(Anno, -1)
    } else {
      stop("the first row should be the description of data")
    }
    if (grep("Download\\s+complete", Anno[length(Anno)])) {
      Anno <- head(Anno, -1)
    } else {
      stop("the last row should be the tag of data")
    }
    tmpf<-tempfile()
    writeLines(Anno,con = tmpf)
    df<-as.data.frame(fread(tmpf,sep='\t',header = FALSE,skip = 1L))
    colnames(df)<-cn
    return(df)
  }
  if (!is.null(login.info)&!is.null(webkey))
      public <- FALSE
    else
      public <- TRUE
    if(is.null(webkey))
      status <- checkMgrastMetagenome(metagenome.id = MetagenomeID, login.info = login.info, public = public)
    else
      status<-TRUE
    if (status){
      MetagenomeID <- paste0("mgm", MetagenomeID)
      server.resource <- "http://api.metagenomics.anl.gov/1/annotation/similarity/"
      server.resource <- paste0(server.resource,MetagenomeID)
      message(paste("Loading the annotations form MG-RAST of", MetagenomeID), domain = NA)
      message("The time spent in this step is proportional to the total amount of remote data...")
      if(!is.null(webkey)){
        param <- list(source = resource["source"], type = resource["type"], evalue = evalue,
                      identity = identity, length = length, auth = webkey)
      } else if (!is.null(login.info)){
        webkey <- login.info["webkey"]
        param <- list(source = resource["source"], type = resource["type"], evalue = evalue,
                      identity = identity, length = length, auth = webkey)
      }else{
        param <- list(source = resource["source"], type = resource["type"], evalue = evalue,
                      identity = identity, length = length)
      }
      anno <- tryCatch(
        getForm(server.resource,
                .params = param,
                .opts=list(noprogress = FALSE,
                           progressfunction = function(down,up){
                             cat(paste('\r', "loading", paste0(round(down[2]/1024^2, 2),"MB"),"..."))
                             flush.console()
                           })
        ),
        error = function(e) {
          msg <- conditionMessage(e)
          structure(msg, class = "try-error")
        }
      )
      if (inherits(anno, "try-error")){
        warning(anno)
        return(FALSE)
      }
      invalid.source <- which(grepl("Invalid\\s+ontology\\s+source", anno))
      if (length(invalid.source))
        stop("invalid ontology source")
      if (length(which(grepl("insufficient\\s+permissions", anno))))
        stop("invalid webkey")
      anno<-checkAnno(anno)
      return(anno)
      cat("\n",MetagenomeID, "annotation data loading completed")
    }else{
      return(NULL)
    }
}

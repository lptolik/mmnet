mapAnnotoKEGGPathway <-
function(KEGGAnno,map.color="blue",browse=FALSE)
{
	url <- sprintf("http://www.kegg.jp/kegg-bin/show_pathway?%s/", 
        	"ec01100")
   	segs <- sprintf("%s%%09%s,%s",KEGGAnno, "", map.color)
        url <- sprintf("%s%s", url, paste(segs, collapse = "/"))
        if (browse)
            browseURL(url)
        return(url)
 }

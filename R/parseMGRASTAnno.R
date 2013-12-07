parseMGRASTAnno <-function(MGRAST.ontology)
{
	#---------------------extracting ec number of enzyme-------------------------------------------------
        message("Getting KEGG ontology annotation,  ...", domain = NA)
        ec.info=gregexpr('(EC:.*)\\s+K\\d+\\s+KO$',MGRAST.ontology,perl=TRUE)
        if (length(ec.info)<1)
            stop("the input MGRAST.ontology need to be the format of functional ontology from the MGRAST") 
        ec.info2=ec.info[ec.info!=-1]
	ec.row=which(ec.info!=-1)
	rast2=MGRAST.ontology[ec.row]
	EC=mapply(function(str,greg)substring(str,attr(greg,'capture.start')[,1],
                   attr(greg,'capture.start')[,1]+attr(greg,'capture.length')[,1]-1),
                   rast2,ec.info2)
	names(EC)=NULL
	ec2=lapply(EC,function(x)gsub('\\]','',x))
	#-----------separate the enzymes mapped to the same read-----------------------
        message("Separatingh the enzymes mapped to the same read,  ...", domain = NA)
	ec2=lapply(ec2,function(x)gsub('\\s',' EC:',x))
	split.ec=mapply(function(x)strsplit(x,'\\s'),ec2)

	unlist.ec=unlist(split.ec)
	#unique.ec=unique(unlist.ec)
	unique.ec=as.data.frame(table(unlist.ec))
        return(as.vector(unique.ec[,1]))
}

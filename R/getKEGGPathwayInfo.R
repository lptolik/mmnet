getKEGGPathwayInfo<- function ()
{
      #-----------------------get substrate from kgml-------------------
      substrateGet <- function(ECPathwayKgml)
     {    
          spkr=strsplit(ECPathwayKgml,'\\s+</reaction>\\n')
          if (length(spkr[[1]])==1)
              stop("the reading KGML file may be fault")
          spkr=spkr[[1]][-length(spkr[[1]])]
          substrate=list()
          substrate.info=mapply(function(x){gregexpr('<substrate\\s+id=\\"\\d+\\"\\s+name=\\"(.*)\\"/>',
                               x,perl=TRUE)},spkr)
          for (i in 1:length(spkr)){
              if (length(attr(substrate.info[[i]],'capture.start')==1)){
              substrate[[i]]=substring(spkr[i],attr(substrate.info[[i]],'capture.start'),
                                       attr(substrate.info[[i]],'capture.start')+
                                       attr(substrate.info[[i]],'capture.length')-1)
     	      }
     	      else{ 
                  for (j in 1:length(attr(substrate.info[[i]],'capture.start'))){
                  substrate[[i]][j]=substring(spkr[i],attr(substrate.info[[i]][j],'capture.start'),
                                             attr(substrate.info[[i]][j],'capture.start')+
                                             attr(substrate.info[[i]][j],'capture.length')-1)
         	  }
     	      }
	 }
         return(substrate)
     }

     #------------------------get product from kgml---------------------------------
      productGet <- function(ECPathwayKgml)
     {
          spkr=strsplit(ECPathwayKgml,'\\s+</reaction>\\n')
          if (length(spkr[[1]])==1)
              stop("the reading KGML file may be fault")
          spkr=spkr[[1]][-length(spkr[[1]])]
          product=list()
          product.info=mapply(function(x){gregexpr('<product\\s+id=\\"\\d+\\"\\s+name=\\"(.*)\\"/>',
                              x,perl=TRUE)},spkr)
          for (i in 1:length(spkr)){
              if (length(attr(product.info[[i]],'capture.start')==1)){
                  product[[i]]=substring(spkr[i],attr(product.info[[i]],'capture.start'),
                                         attr(product.info[[i]],'capture.start')+
                                         attr(product.info[[i]],'capture.length')-1)
              }
     	      else{   
                  for (j in 1:length(attr(product.info[[i]],'capture.start'))){
                      product[[i]][j]=substring(spkr[i],attr(product.info[[i]][j],'capture.start'),
                                                attr(product.info[[i]][j],'capture.start')+
                                                attr(product.info[[i]][j],'capture.length')-1)
                  }
     	     }
          }
          return(product)
      }

      ECMetabolites <- function(ec,substrate,product,rn.type)
      {
          if (!length(ec)==length(substrate) && length(product)==length(substrate) && length(product)!=length(rn.type))
              stop(" the length of ec, substrate, product and rn.type must be the same") 
          paste.metabolites=function(x){
              b=x[1] 
              if(length(x)>1){
     	         for ( i in 2:length(x))
                     b=paste(b,x[i])
 	      }
          return(b)
          }

          #-------------a reaction may correspond to several ecs------------------------
          ec=unlist(ec)
          rn.type=unlist(rn.type)
          paste.substrate=mapply(function(x)mapply(paste.metabolites,x),substrate)
          paste.product=mapply(function(x)mapply(paste.metabolites,x),product)
          unlist.substrate=unlist(paste.substrate)
          unlist.product=unlist(paste.product)
          Multiecinfo=gregexpr('(\\s)',ec,perl=TRUE)
          tmp=mapply(function(x)isTRUE(x[1]!=-1),Multiecinfo)
          Multiec.index=which(tmp==TRUE)
          Multiec=strsplit(ec[Multiec.index],'\\s')
          tmpec=ec
          tmpec[Multiec.index]=mapply(function(x)x[1],Multiec)
          Multiec2=(mapply(function(x)x[-1],Multiec))
          Multiec3=unlist(Multiec2)
          tmpec=append(tmpec,Multiec3) 
          Multimetabolites.index=rep(Multiec.index,mapply(length,Multiec2))
          Multisubstrate=unlist.substrate[Multimetabolites.index]
          Multiproduct=unlist.product[Multimetabolites.index]
          Multi.rn.type=rn.type[Multimetabolites.index]
          append.substrate=append(unlist.substrate,Multisubstrate)
          append.product=append(unlist.product,Multiproduct)
          append.rn.type=append(rn.type,Multi.rn.type)

          #-------------------processing of reversible reaction -----------------------
          message("Processing of reversible reaction  ...", domain = NA)
          reverse.index=which(append.rn.type=="reversible")
          reverse.substrate=append.substrate
          reverse.substrate[reverse.index]=paste(append.substrate[reverse.index],append.product[reverse.index])
          reverse.product=append.product
          reverse.product[reverse.index]=paste(append.substrate[reverse.index],append.product[reverse.index])

          #----------------- the uique of ec-------------------------------------------
          message("Merging the substrates and products of the same enzyme  ...", domain = NA)
          delete.no=which(duplicated(tmpec)==TRUE)
          dup=unique(tmpec[which(isUnique(tmpec)==FALSE)])
          dup.no=sapply(dup,function(x)which(tmpec==x))

          #------merge the substrates of the same ec------------------
          merging.substrate=reverse.substrate
          for (i in 1:length(dup.no)){
              item=unlist(dup.no[[i]])
              j=2      
              while(j<=length(item)){
                  merging.substrate[item[1]]=paste(merging.substrate[item[1]],merging.substrate[item[j]])
           	  j=j+1	 
              }
          }
          item=c(1:length(merging.substrate))
          item=item[-delete.no]
          substrate=mapply(function(x)merging.substrate[x],item)

          #-------------merge the products of the same ec---------------
          merging.product=reverse.product
          for (i in 1:length(dup.no)){
              item=unlist(dup.no[[i]])
              j=2      
              while(j<=length(item)){
            	  merging.product[item[1]]=paste(merging.product[item[1]],merging.product[item[j]])
           	  j=j+1	 
              }
          }
          item=c(1:length(merging.product))
          item=item[-delete.no]
          product=mapply(function(x)merging.product[x],item)
     	  ec=tmpec[-delete.no]

          #------delete the same metabolites for the each ec-----------------
          message("Deleting the same metabolites for each enzyme  ...", domain = NA)
          dlt.dup.metabolites=function(x){
              str.metabio=strsplit(x,'\\s')
              y=unique(unlist(str.metabio))
          }
          unique.substrate=mapply(dlt.dup.metabolites,substrate)
          paste.substrate2=lapply(unique.substrate,paste.metabolites)
          substrate=unlist(paste.substrate2,use.names=F)
          unique.product=mapply(dlt.dup.metabolites,product)
          paste.product2=lapply(unique.product,paste.metabolites)
          product=unlist(paste.product2,use.names=F)
 	  return(list(ec=ec,substrate=substrate,product=product))
      }
      
      #-------download the KGML file of KEGG ec reference pathway except the global metabolic pathway-------
      pathway.ec=keggList("pathway/ec")
      ec.info=regexpr('ec\\d+',names(pathway.ec),perl=TRUE)
      ec=noquote(substring(names(pathway.ec),ec.info,ec.info+attr(ec.info,'match.length')-1))
      global.ec=ec %in% c('ec01100','ec01110','ec01120')
      new.ec=ec[!global.ec]
      kgml=""
      message("Getting the KGML files of the KEGG ec reference pathway with the KEGG API  ...", domain = NA)
      for (ec.no in new.ec){
          temp=try(keggGet(ec.no,"kgml"),silent=T)
          if (class(temp)!="try-error"){
              kgml=paste(kgml,temp)
          }
      }     
      message("Parsing the KGML files to get the information of enzyme and metabolites  ...", domain = NA)
      spkr=strsplit(kgml,'<pathway')
      spkr=spkr[[1]][-1]
      pathway.no.info=lapply(spkr,function(x){gregexpr('name=\\"path:(ec\\d+)\\"\\s+org=',
			     x,perl=TRUE)})
      pathway.no=list()
      for (i in 1:length(pathway.no.info)){
          if (attr(pathway.no.info[[i]][[1]],"match.length")>0){
	      pathway.no[i]=substring(spkr[i],attr(pathway.no.info[[i]][[1]],
				     'capture.start'),attr(pathway.no.info[[i]][[1]],
				     'capture.start')+attr(pathway.no.info[[i]][[1]],
				    'capture.length')-1)  
	  }
      }

      #-------------------------------------info of enzyme, reaction and metabolites------------------------
      message("Extracting the information of enzyme, reaction  ...", domain = NA)
      spkr=as.list(spkr)
      ec.info=sapply(spkr,function(x)gregexpr('id=\\"(\\d+)\\"\\s+name=\\"(ec:.*)\\"\\s+type',x,perl=TRUE))
      ec.id=mapply(function(str,greg)substring(str,attr(greg,'capture.start')[,1],
                   attr(greg,'capture.start')[,1]+attr(greg,'capture.length')[,1]-1),
                   spkr,ec.info)
      ec.no=mapply(function(str,greg)substring(str,attr(greg,'capture.start')[,2],
                   attr(greg,'capture.start')[,2]+attr(greg,'capture.length')[,2]-1),
                   spkr,ec.info)
      rn.info=sapply(spkr,function(x)gregexpr('id=\\"(\\d+)\\"\\s+name=\\"(rn:.*)\\"\\s+type=\\"(\\w+)\\">',
                       x,perl=TRUE))
      rn.id=mapply(function(str,greg)substring(str,attr(greg,'capture.start')[,1],
                   attr(greg,'capture.start')[,1]+attr(greg,'capture.length')[,1]-1),
                   spkr,rn.info)
      rn.type=mapply(function(str,greg)substring(str,attr(greg,'capture.start')[,3],
                     attr(greg,'capture.start')[,3]+attr(greg,'capture.length')[,3]-1),
                     spkr,rn.info) 
      ec.id2=mapply(intersect,ec.id,rn.id)     
      #--------------------------make enzyme and reaction correspondence ----------------------------
      message("Making enzyme, reaction and metabolites correspondence  ...", domain = NA)
      select.log=mapply(is.element,ec.id,ec.id2)
      select.index=mapply(function(x)which(x==TRUE),select.log)
      delete.index=which(mapply(length,ec.id2)==0)
      rn.id[delete.index]=NULL
      ec.id2[delete.index]=NULL
      spkr[delete.index]=NULL
      select.index[delete.index]=NULL
      ec.no[delete.index]=NULL
      rn.type[delete.index]=NULL
      ec.no2=mapply(function(a,b)a[b],ec.no,select.index)
      fsort.id=function(a,b)mapply(function(x)which(x==a),b)
      sort.id=mapply(fsort.id,rn.id,ec.id2)
      message("Extracting the substrates and products ...", domain = NA)  
      substrate=mapply(substrateGet,spkr)
      product=mapply(productGet,spkr)
      substrate.sort=mapply(function(substrate,sort.id)substrate[sort.id],substrate,sort.id)
      product.sort=mapply(function(product,sort.id)product[sort.id],product,sort.id)
      rn.type.sort=mapply(function(rn.type,sort.id)rn.type[sort.id],rn.type,sort.id)
      ec=ec.no2
      substrate=substrate.sort
      product=product.sort
      rn.type=rn.type.sort
      KEGGPathwayInfo=ECMetabolites(ec,substrate,product,rn.type)
      message("Done ...", domain = NA)
      #return(KEGGPathwayInfo) 
      path=system.file("data",package="mmnet")
      save(KEGGPathwayInfo,file=sprintf("%s/%s",path,"KEGGPathwayInfo.rda"))
}
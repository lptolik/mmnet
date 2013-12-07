constructMetabolicNetwork <-
function(KEGGPathwayInfo,MGRAST.KEGGAnno,degree.cutoff=0)
{
      if (length(unique(mapply(length,KEGGPathwayInfo)))!=1)
          stop("the input KEGGPathwayInfo is wrong") 
      ec=KEGGPathwayInfo$ec
      product=KEGGPathwayInfo$product
      substrate=KEGGPathwayInfo$substrate
      ec=gsub('ec','EC',ec)
      message("Constructing nodes of the network  ...", domain = NA)
      node1=intersect(ec,MGRAST.KEGGAnno)
      node1.index=mapply(function(node1)which(ec==node1),node1)
      sub.substrate=substrate[node1.index]
      sub.product=product[node1.index]
      prepro.metabo=function(x){
                 str.metabo=strsplit(x,'\\s')
                 str.metabo=unlist(str.metabo)
      }
      str.product=lapply(sub.product,prepro.metabo)
      str.substrate=lapply(sub.substrate,prepro.metabo)
      message("Constructing edges of the network  ...", domain = NA)
      edge.matrix=mapply(function(y)mapply(function(x)intersect(y,x),str.substrate),str.product)      
      Getedge=function(x){length(x[[1]])}
      edge=apply(edge.matrix,c(1,2),Getedge)
      diag(edge)=0
      row.index=apply(edge,2,function(x)which(x!=0))
      edge[edge>1]=1
      edge=t(edge)
      NW=graph.adjacency(edge)
      bad.vs<-V(NW)[degree(NW)<=degree.cutoff]
      NW<-delete.vertices(NW, bad.vs)
}

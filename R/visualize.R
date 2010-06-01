##############################################################
## plot Graph
plotAnn<-function(pathway,graphList,ann,gotoKEGG=FALSE){
     require(Rgraphviz)
   
     nAttrs<-list()
      geneList<-ann[[pathway]]$annGeneList
      ecList<-getEnzymeFromGene(geneList)
      z<-rep("red",length(ecList))
      names(z)<-ecList
      nAttrs$color<-z
      plot(graphList[[pathway]],"neato",nodeAttrs=nAttrs)
      if(gotoKEGG==TRUE){
      kgeneList<-getKidFromOid(geneList)
      org<-getOrgAndIdType()[1]
      pathwayId<-substring(pathway,6,10)  
	  temp<- paste ( c ( paste ( org , pathwayId , sep = "" ) , kgeneList ) , sep = "" , collapse = "+" )
      url <- paste("http://www.genome.ad.jp/dbget-bin/show_pathway?", temp , sep = "")
      browseURL(url)
      }
}
##############################################################
##new! plot Graph
plotKOAnn<-function(pathway,graphList,ann,gotoKEGG=FALSE){
     require(Rgraphviz)
   
     nAttrs<-list()
      geneList<-ann[[pathway]]$annGeneList
      KOList<-getKOFromGene(geneList)
      z<-rep("red",length(KOList))
      names(z)<-KOList
      nAttrs$color<-z
      plot(graphList[[pathway]],"neato",nodeAttrs=nAttrs)
      if(gotoKEGG==TRUE){
      kgeneList<-getKidFromOid(geneList)
      org<-getOrgAndIdType()[1]
      pathwayId<-substring(pathway,6,10)  
	  temp<- paste ( c ( paste ( org , pathwayId , sep = "" ) , kgeneList ) , sep = "" , collapse = "+" )
      url <- paste("http://www.genome.ad.jp/dbget-bin/show_pathway?", temp , sep = "")
      browseURL(url)
      }
}
############################################################
##go to html
gotoKEGG<-function(pathway,ann){
      geneList<-ann[[pathway]]$annGeneList
      kgeneList<-getKidFromOid(geneList)
      org<-getOrgAndIdType()[1]
      pathwayId<-substring(pathway,6,10)  
	  temp<- paste ( c ( paste ( org , pathwayId , sep = "" ) , kgeneList ) , sep = "" , collapse = "+" )
      url <- paste("http://www.genome.ad.jp/dbget-bin/show_pathway?", temp , sep = "")
      browseURL(url)
}

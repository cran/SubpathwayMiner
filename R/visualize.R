##############################################################
## plot Graph
plotAnn<-function(pathway,graphList,ann,gotoKEGG=FALSE){
     library(Rgraphviz)
   
     nAttrs<-list()
      geneList<-ann[[pathway]]$annGeneList
      ecList<-getEnzymeFromGene(geneList)
      z<-rep("red",length(ecList))
      names(z)<-ecList
      nAttrs$color<-z
      plot(graphList[[pathway]],"neato",nodeAttrs=nAttrs)
      if(gotoKEGG==TRUE){
           kgeneList<-getKidFromOid(geneList)
           s<-character()
           for(i in 1:length(kgeneList)){
                  s<-paste(s,kgeneList[i],sep="+")
           }
           org<-getOrgAndIdType()[1]
           pathwayId<-substring(pathway,6,10)
           url<-paste("http://www.genome.ad.jp/dbget-bin/show_pathway?",org,pathwayId,s,sep="")
           browseURL(url)
      }
}
##############################################################
##new! plot Graph
plotKOAnn<-function(pathway,graphList,ann,gotoKEGG=FALSE){
     library(Rgraphviz)
   
     nAttrs<-list()
      geneList<-ann[[pathway]]$annGeneList
      KOList<-getKOFromGene(geneList)
      z<-rep("red",length(KOList))
      names(z)<-KOList
      nAttrs$color<-z
      plot(graphList[[pathway]],"neato",nodeAttrs=nAttrs)
      if(gotoKEGG==TRUE){
           kgeneList<-getKidFromOid(geneList)
           s<-character()
           for(i in 1:length(kgeneList)){
                  s<-paste(s,kgeneList[i],sep="+")
           }
           org<-getOrgAndIdType()[1]
           pathwayId<-substring(pathway,6,10)
           url<-paste("http://www.genome.ad.jp/dbget-bin/show_pathway?",org,pathwayId,s,sep="")
           browseURL(url)
      }
}
############################################################
##go to html
gotoKEGG<-function(pathway,ann){
      library(Rgraphviz)
      geneList<-ann[[pathway]]$annGeneList
      kgeneList<-getKidFromOid(geneList)
      s<-character()
      for(i in 1:length(kgeneList)){
            s<-paste(s,kgeneList[i],sep="+")
      }
      org<-getOrgAndIdType()[1]
      pathwayId<-substring(pathway,6,10)
      url<-paste("http://www.genome.ad.jp/dbget-bin/show_pathway?",org,pathwayId,s,sep="")
      browseURL(url)
}

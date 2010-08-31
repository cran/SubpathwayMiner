 
####################################################################
##get annotation
getAnn<-function(geneList,background=getDefaultBackground(),
   order="pvalue",decreasing=FALSE,graphList=getDefaultGraph()){
      if(typeof(geneList)!="character"){
	  print("warning: your geneList must be 'character' vector. Because the type of your current geneList is not correct, it has been conveted arbitrarily using the function as.character().")
	  geneList<-as.character(geneList)
	  }
      if(!exists("ke2g")) initialize_ke2g()
      graphList<-graphList[sapply(graphList,function(x) length(x)>0)]     
      keggpathid2name<-get("keggpathid2name",envir=ke2g)

      annList<-list()
      for(i in 1:length(graphList)){
            ann<-list(pathwayName="not known",annGeneList=character(),annBgGeneList=character(),annGeneNumber=0,
                      annBgNumber=0,geneNumber=0,bgNumber=0,pvalue=1,qvalue=1)
            if(class(graphList[[i]])=="character"){
                  graphGeneList<-getGeneFromPathway(graphList[[i]])
				  pathwayName<-getPNameFromPId(graphList[[i]])
            }
            else if(class(graphList[[i]])=="graphNEL"){
                  graphGeneList<-getGeneFromEnzyme(nodes(graphList[[i]]))   
                  pathwayName<-keggpathid2name[[substring(names(graphList)[i],6,10)]]				  
            }            
            annotatedGeneList<-intersect(graphGeneList,geneList)
            annotatedBackgroundList<-intersect(graphGeneList,background)

            #pathwayName<-keggpathid2name[names(graphList)[i]]
            if(length(pathwayName)!=0)
                ann$pathwayName<-pathwayName
            ann$annGeneList<-annotatedGeneList   
            ann$annBgGeneList<-annotatedBackgroundList
         
            ann$annGeneNumber<-length(annotatedGeneList)
            ann$annBgNumber<-length(annotatedBackgroundList)

            ann$geneNumber<-length(geneList)
            ann$bgNumber<-length(background)

            ann$pvalue<-1-phyper(ann$annGeneNumber-1,ann$annBgNumber,
                 ann$bgNumber-ann$annBgNumber,ann$geneNumber)
            
            annList[[i]]<-ann
      } 
      qvalueList<-fdrtool(sapply(annList,function(x) x$pvalue), 
                   statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
      for(i in 1:length(annList)){
            annList[[i]]$qvalue<-qvalueList[i]
      }
      names(annList)<-names(graphList)
      annList<-annList[sapply(annList,function(x) x$annGeneNumber>0)]
      annList<-annList[order(sapply(annList,function(x) x[[order]]),decreasing=decreasing)]

      return(annList)
}
####################################################################
##get metabolic pathway annotation
getMpAnn<-function(geneList,background=getDefaultBackground(),
   order="pvalue",decreasing=FALSE){
getAnn(geneList=geneList,background=background,
   order=order,decreasing=decreasing,graphList=getDefaultUndirectedGraph())

}
####################################################################
##get sub-pathway annotation of metabolic pathways 
getKcsmpAnn<-function(geneList,background=getDefaultBackground(),k=4,
   order="pvalue",decreasing=FALSE){

      subGraph<-getKcSubGraph(k,graphList=getDefaultUndirectedGraph())
      ann<-getAnn(geneList=geneList,background=background,order=order,decreasing=decreasing,graphList=subGraph)
      return(ann)
}
####################################################################
##cutoff annotation
cutoffAnn<-function(ann,type="pvalue",operate="<=",cutoff=0.01){
      if(operate=="<"){
            ann<-ann[sapply(ann,function(x) x[[type]])<cutoff]
      }
      else if(operate==">"){
            ann<-ann[sapply(ann,function(x) x[[type]])>cutoff]
      }
      else if(operate=="<="){
            ann<-ann[sapply(ann,function(x) x[[type]])<=cutoff]
      }
      else if(operate==">="){
            ann<-ann[sapply(ann,function(x) x[[type]])>=cutoff]
      }
      return(ann)
}

#####################################################################
##print Ann
printAnn<-function(ann,detail=FALSE){
	  if(detail==FALSE){
      pathwayName<-sapply(ann,function(x) x$pathwayName)
      annGeneRatio<-sapply(ann,function(x) paste(x$annGeneNumber,x$geneNumber,sep="/"))
      annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
      pvalue<-sapply(ann,function(x) x$pvalue)
      qvalue<-sapply(ann,function(x) x$qvalue)
      ann.data.frame<-as.data.frame(cbind(pathwayName,annGeneRatio,
                             annBgRatio,pvalue,qvalue))
	  }
	  else{	  
	  pathwayName<-sapply(ann,function(x) x$pathwayName)
	  annGeneList<-sapply(ann, function(x){ paste(x$annBgGeneList,collapse=";") })
      annBgGeneList<-sapply(ann, function(x){ paste(x$annBgGeneList,collapse=";")})
	  annGeneRatio<-sapply(ann,function(x) paste(x$annGeneNumber,x$geneNumber,sep="/"))
      annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
      pvalue<-sapply(ann,function(x) x$pvalue)
      qvalue<-sapply(ann,function(x) x$qvalue)
      ann.data.frame<-as.data.frame(cbind(pathwayName,annGeneRatio,
                             annBgRatio,pvalue,qvalue,annGeneList,annBgGeneList))
							 
	  }
      return(ann.data.frame)
}


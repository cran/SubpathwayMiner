
#############################################################
##get kegg gene list from gene list
getKGeneFromGene<-function(geneList){
	  geneList<-as.character(geneList)
      if(!exists("ke2g")) initialize_ke2g()
	  keggGene2gene<-get("keggGene2gene",envir=ke2g)
keggGeneList<-unique(as.character(keggGene2gene[as.character(keggGene2gene[,2]) %in% paste(getOrgAndIdType()[2],geneList,sep=":"),1]))
      return(keggGeneList)
}
#keggGeneList<-getKeggGeneFromGene(geneList)
#############################################################
##get gene list from kegg gene list
getGeneFromKGene<-function(keggGeneList){
	  keggGeneList<-as.character(keggGeneList)
      if(!exists("ke2g")) initialize_ke2g()
      keggGene2gene<-get("keggGene2gene",envir=ke2g)
      geneList<-unique(as.character(sapply(strsplit(as.character(keggGene2gene[as.character(keggGene2gene[,1]) %in% keggGeneList,2]),":"),function(x) return (x[2]))))
      return(geneList)
}
#getGeneFromKeggGene(keggGeneList)
#############################################################
##get enzyme list from gene List
getEnzymeFromKGene<-function(keggGeneList){
	  keggGeneList<-as.character(keggGeneList)
      if(!exists("ke2g")) initialize_ke2g()
	  gene2ec<-get("gene2ec",envir=ke2g)
      enzymeList<-unique(as.character(gene2ec[as.character(gene2ec[,1]) %in% keggGeneList,2]))
      return(enzymeList)
}
#enzymeList<-getEnzymeFromKeggGene(keggGeneList)
#############################################################
##get gene list from enzyme List
getKGeneFromEnzyme<-function(enzymeList){
	  enzymeList<-as.character(enzymeList)
      if(!exists("ke2g")) initialize_ke2g()
	  gene2ec<-get("gene2ec",envir=ke2g)
      keggGeneList<-unique(as.character(gene2ec[as.character(gene2ec[,2]) %in% enzymeList,1]))
      return(keggGeneList)
}
#getKeggGeneFromEnzyme(enzymeList)
#############################################################
##new! get KO list from gene List
getKOFromKGene<-function(keggGeneList){
	  keggGeneList<-as.character(keggGeneList)
      if(!exists("ke2g")) initialize_ke2g()
	  gene2ko<-get("gene2ko",envir=ke2g)
      KOList<-unique(as.character(gene2ko[as.character(gene2ko[,1]) %in% keggGeneList,2]))
      return(KOList)
}
#KOList<-getKOFromKeggGene(keggGeneList)
#############################################################
##new! get gene list from KO List
getKGeneFromKO<-function(KOList){
	  KOList<-as.character(KOList)
      if(!exists("ke2g")) initialize_ke2g()
	  gene2ko<-get("gene2ko",envir=ke2g)
      keggGeneList<-unique(as.character(gene2ko[as.character(gene2ko[,2]) %in% KOList,1]))
      return(keggGeneList)
}
#getKeggGeneFromKO(KOList)
#############################################################
##new! get pathway from gene List
getPathwayFromKGene<-function(keggGeneList){
	  keggGeneList<-as.character(keggGeneList)
      if(!exists("ke2g")) initialize_ke2g()
	  gene2path<-get("gene2path",envir=ke2g)
      pathwayList<-unique(as.character(gene2path[as.character(gene2path[,1]) %in% keggGeneList,2]))
	  pathwayList<-paste("path:",substring(pathwayList,9),sep="")
      return(pathwayList)
}
#KOList<-getKOFromKeggGene(keggGeneList)
#############################################################
##new! get gene list from pathway List
getKGeneFromPathway<-function(pathwayList){
	  pathwayList<-as.character(pathwayList)
      if(!exists("ke2g")) initialize_ke2g()
	  gene2path<-get("gene2path",envir=ke2g)  
      keggGeneList<-unique(as.character(gene2path[as.character(gene2path[,2]) %in% paste("path:",getOrgAndIdType()[1],substring(pathwayList,6),sep=""),1]))
      return(keggGeneList)
}

#############################################################
getPNameFromPId<-function(PIdList){
	  PIdList<-as.character(PIdList)
      if(!exists("ke2g")) initialize_ke2g()
	  map_title<-get("map_title",envir=ke2g)
      PNameList<-unique(as.character(map_title[as.character(map_title[,1]) %in% substring(PIdList,6),2]))
      return(PNameList)
}
#############################################################
getPIdFromPName<-function(PNameList){
	  PNameList<-as.character(PNameList)
      if(!exists("ke2g")) initialize_ke2g()
	  map_title<-get("map_title",envir=ke2g)
      PIdList<-unique(as.character(map_title[as.character(map_title[,2]) %in% PNameList,1]))
      return(PIdList)
}
#####################################################################
getEnzymeFromGene<-function(geneList){
      return(getEnzymeFromKGene(getKGeneFromGene(geneList)))
}
getGeneFromEnzyme<-function(enzymeList){
      return(getGeneFromKGene(getKGeneFromEnzyme(enzymeList)))
}
getKOFromGene<-function(geneList){
      return(getKOFromKGene(getKGeneFromGene(geneList)))
}
getGeneFromKO<-function(KOList){
      return(getGeneFromKGene(getKGeneFromKO(KOList)))
}
getPathwayFromGene<-function(geneList){
      return(getPathwayFromKGene(getKGeneFromGene(geneList)))
}
getGeneFromPathway<-function(pathwayList){
      return(getGeneFromKGene(getKGeneFromPathway(pathwayList)))
}
#############################################################
##get defaultMetabolic pathway
getDefaultMetabolicPathway<-function(){
      if(!exists("ke2g")) initialize_ke2g()
      mpidList<-get("mpidList",envir=ke2g)
      return(mpidList)
}
#############################################################
##new! get default KO pathway
getDefaultKOPathway<-function(){
      if(!exists("ke2g")) initialize_ke2g()
      kpidList<-get("kpidList",envir=ke2g)
      #kpidList<-get("keggpathid2name",envir=ke2g)
      #kpidList<-paste("path:",names(kpidList),sep="")
      return(kpidList)
}

####################################################################
##
getDefaultBackground<-function(){
      if(!exists("ke2g")) initialize_ke2g()
	  keggGene2gene<-get("keggGene2gene",envir=ke2g) 
	  background<-unique(as.character(keggGene2gene[,2]))
	  newBackground<-sapply(strsplit(background,":"),function(x) x[2])
      return(newBackground)
}
#####################################################################
##
getDefaultGraph<-function(){
      if(!exists("ke2g")) initialize_ke2g()
      gene2path<-get("gene2path",envir=ke2g)
	  graphList<-unique(paste("path:",substring(gene2path[,2],9),sep=""))
	  names(graphList)<-graphList
      return(graphList)   
}
#####################################################################
##
getDefaultUndirectedGraph<-function(){
      if(!exists("ke2g")) initialize_ke2g()
      uGraph<-get("uGraph",envir=ke2g)
      return(uGraph)       
}
#####################################################################
##new!
getDefaultKOUndirectedGraph<-function(){
      if(!exists("ke2g")) initialize_ke2g()
      KOuGraph<-get("KOuGraph",envir=ke2g)
      return(KOuGraph)       
}
#####################################################################
getOrgAndIdType<-function(){
      if(!exists("ke2g")) initialize_ke2g()
      orgAndIdType<-get("orgAndIdType",envir=ke2g)
      return(orgAndIdType)    
}
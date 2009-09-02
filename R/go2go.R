#############################################################
##get pathway name from id
getPathwayNameFromId<-function(pid){
      if(!exists("ke2g")) initialize()
      pid<-substring(pid,6,10)
      keggpathid2name<-get("keggpathid2name",envir=ke2g)
      name<-unlist(keggpathid2name[pid])
      return(name)
}
#############################################################
##get kid from oid
getKidFromOid<-function(oid){
      if(!exists("ke2g")) initialize()
      oid2kid<-get("oid2kid",envir=ke2g)
      kid<-unique(unlist(oid2kid[oid]))
      return(kid)
}
#############################################################
##get oid list from kid
getOidFromKid<-function(kid){
      if(!exists("ke2g")) initialize()
      kid2oid<-get("kid2oid",envir=ke2g)
      oid<-unique(unlist(kid2oid[kid]))
      return(oid)
}
#############################################################
##get enzyme list from gene List
getEnzymeFromGene<-function(geneList){
      if(!exists("ke2g")) initialize()
      gene2ec<-get("gene2ec",envir=ke2g)
      enzymeList<-unique(unlist(gene2ec[geneList]))
      return(enzymeList)
}
#############################################################
##get gene list from enzyme List
getGeneFromEnzyme<-function(enzymeList){
      if(!exists("ke2g")) initialize()
      ec2gene<-get("ec2gene",envir=ke2g)
      geneList<-unique(unlist(ec2gene[enzymeList]))
      return(geneList)
}
#############################################################
##new! get KO list from gene List
getKOFromGene<-function(geneList){
      if(!exists("ke2g")) initialize()
      gene2KO<-get("gene2KO",envir=ke2g)
      KOList<-unique(unlist(gene2KO[geneList]))
      return(KOList)
}
#############################################################
##new! get gene list from KO List
getGeneFromKO<-function(KOList){
      if(!exists("ke2g")) initialize()
      KO2gene<-get("KO2gene",envir=ke2g)
      geneList<-unique(unlist(KO2gene[KOList]))
      return(geneList)
}
#############################################################
##get pathway list from gene list
getPathwayFromGene<-function(geneList){
      if(!exists("ke2g")) initialize()
      gene2path<-get("gene2path",envir=ke2g)
      pathwayList<-unique(unlist(gene2path[geneList]))
      return(pathwayList)
}
#############################################################
##get gene list from pathway list
getGeneFromPathway<-function(pathwayList){
      if(!exists("ke2g")) initialize()
      path2gene<-get("path2gene",envir=ke2g)
      geneList<-unique(unlist(path2gene[pathwayList]))
      return(geneList)
}
#############################################################
##get defaultMetabolic pathway
getDefaultMetabolicPathway<-function(){
      if(!exists("ke2g")) initialize()
      mpidList<-get("mpidList",envir=ke2g)
      return(mpidList)
}
#############################################################
##new! get default KO pathway
getDefaultKOPathway<-function(){
      if(!exists("ke2g")) initialize()
      kpidList<-get("kpidList",envir=ke2g)
      #kpidList<-get("keggpathid2name",envir=ke2g)
      #kpidList<-paste("path:",names(kpidList),sep="")
      return(kpidList)
}

####################################################################
##
getDefaultBackground<-function(){
      if(!exists("ke2g")) initialize()
      background<-get("background",envir=ke2g)
      return(background)
}
#####################################################################
##
getDefaultGraph<-function(){
      if(!exists("ke2g")) initialize()
      path2gene<-get("path2gene",envir=ke2g)
      return(path2gene)   
}
#####################################################################
##
getDefaultUndirectedGraph<-function(){
      if(!exists("ke2g")) initialize()
      uGraph<-get("uGraph",envir=ke2g)
      return(uGraph)       
}
#####################################################################
##new!
getDefaultKOUndirectedGraph<-function(){
      if(!exists("ke2g")) initialize()
      KOuGraph<-get("KOuGraph",envir=ke2g)
      return(KOuGraph)       
}
#####################################################################
getOrgAndIdType<-function(){
      if(!exists("ke2g")) initialize()
      orgAndIdType<-get("orgAndIdType",envir=ke2g)
      return(orgAndIdType)    
}
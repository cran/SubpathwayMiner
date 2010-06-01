#############################################################
##get pathway name from id
getPathwayNameFromId<-function(pid){
      if(!exists("ke2g")) initialize_ke2g()
      pid<-substring(pid,6,10)
      keggpathid2name<-get("keggpathid2name",envir=ke2g)
      name<-unlist(keggpathid2name[pid])
      return(name)
}
#############################################################
##get kid from oid
getKidFromOid<-function(oid){
      if(typeof(oid)!="character"){
	  print("warning: your oid must be 'character' vector. Because the type of your current oid is not correct, it has been conveted arbitrarily using the function as.character().")
	  as.character(oid)
	  }
      if(!exists("ke2g")) initialize_ke2g()
      oid2kid<-get("oid2kid",envir=ke2g)
      kid<-unique(unlist(oid2kid[oid]))
      return(kid)
}
#############################################################
##get oid list from kid
getOidFromKid<-function(kid){
      if(typeof(kid)!="character"){
	  print("warning: your kid must be 'character' vector. Because the type of your current kid is not correct, it has been conveted arbitrarily using the function as.character().")
	  as.character(kid)
	  }
      if(!exists("ke2g")) initialize_ke2g()
      kid2oid<-get("kid2oid",envir=ke2g)
      oid<-unique(unlist(kid2oid[kid]))
      return(oid)
}
#############################################################
##get enzyme list from gene List
getEnzymeFromGene<-function(geneList){
      if(!exists("ke2g")) initialize_ke2g()
      gene2ec<-get("gene2ec",envir=ke2g)
      enzymeList<-unique(unlist(gene2ec[geneList]))
      return(enzymeList)
}
#############################################################
##get gene list from enzyme List
getGeneFromEnzyme<-function(enzymeList){
      if(!exists("ke2g")) initialize_ke2g()
      ec2gene<-get("ec2gene",envir=ke2g)
      geneList<-unique(unlist(ec2gene[enzymeList]))
      return(geneList)
}
#############################################################
##new! get KO list from gene List
getKOFromGene<-function(geneList){
      if(!exists("ke2g")) initialize_ke2g()
      gene2KO<-get("gene2KO",envir=ke2g)
      KOList<-unique(unlist(gene2KO[geneList]))
      return(KOList)
}
#############################################################
##new! get gene list from KO List
getGeneFromKO<-function(KOList){
      if(!exists("ke2g")) initialize_ke2g()
      KO2gene<-get("KO2gene",envir=ke2g)
      geneList<-unique(unlist(KO2gene[KOList]))
      return(geneList)
}
#############################################################
##get pathway list from gene list
getPathwayFromGene<-function(geneList){
      if(!exists("ke2g")) initialize_ke2g()
      gene2path<-get("gene2path",envir=ke2g)
      pathwayList<-unique(unlist(gene2path[geneList]))
      return(pathwayList)
}
#############################################################
##get gene list from pathway list
getGeneFromPathway<-function(pathwayList){
      if(!exists("ke2g")) initialize_ke2g()
      path2gene<-get("path2gene",envir=ke2g)
      geneList<-unique(unlist(path2gene[pathwayList]))
      return(geneList)
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
      background<-get("background",envir=ke2g)
      return(background)
}
#####################################################################
##
getDefaultGraph<-function(){
      if(!exists("ke2g")) initialize_ke2g()
      path2gene<-get("path2gene",envir=ke2g)
      return(path2gene)   
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
.first.lib <- function(lib, pkgname, where){
  initialize()   
}
########################################################################
##initialize data
initialize<-function(){
      assign("ke2g",new.env(parent=globalenv()),envir=.GlobalEnv)
      data("hsa_ncbi-geneid")
      assign("keggpathid2name",as.list(KEGGPATHID2NAME),envir=ke2g)
      #assign("orgAndIdType",c("hsa","ncbi-geneid"),envir=ke2g)
      #assign("keggpathid2name",as.list(KEGGPATHID2NAME),envir=ke2g)
      #data(keggpathid2name,envir=ke2g)
      #data(mpidList,envir=ke2g)
      #data(uGraph,envir=ke2g)
      #data(dGraph,envir=ke2g)
      #data(gene2ec,envir=ke2g)  
      #data(ec2gene,envir=ke2g) 
      #data(gene2path,envir=ke2g)
      #data(path2gene,envir=ke2g)
      #data(background,envir=ke2g)
      #data(kid2oid,envir=ke2g)
      #data(oid2kid,envir=ke2g)
}

##########################################################################
##update setting
updateOrgAndIdType<-function(org="hsa",idType="ncbi-geneid",
  path="ftp://ftp.genome.jp/pub/kegg/genes/organisms",verbose=TRUE){
       if(!exists("ke2g")) initialize()
     if(verbose==TRUE){
      print("Note that the programming is time consumming!!!!!!!!!!!!!!")  
      print("download and treat  cross reference identifiers....................")
      }
      kidoid<-getkidoid(org,idType,path)
      assign("kid2oid",kidoid[[1]],envir=ke2g)
       assign("oid2kid",kidoid[[2]],envir=ke2g) 
      if(verbose==TRUE)
      print("download and treat relation between gene and enzyme................")
      geneEnzyme<-getMerge(org,"enzyme",idType,path)
      assign("ec2gene",geneEnzyme[[1]],envir=ke2g)
      assign("gene2ec",geneEnzyme[[2]],envir=ke2g)
      if(verbose==TRUE)
      print("download and treat relation between gene and pathway...............")
      genePathway<-getMerge1(org,"pathway",idType,path)
      assign("path2gene",genePathway[[1]],envir=ke2g)
      assign("gene2path",genePathway[[2]],envir=ke2g)
      assign("background",buildBgGeneList(org,idType,path),envir=ke2g)
      assign("orgAndIdType",c(org,idType),envir=ke2g)
      #assign("keggpathid2name",getPathName(),envir=ke2g)
}
#updateOrgAndIdType("hsa","ncbi-geneid")
###########################################################################
##updateGraph
updateGraphs<-function(pathwayList=getDefaultMetabolicPathway(),
         path="ftp://ftp.genome.jp/pub/kegg/xml/map/",verbose=TRUE){
      library(XML)
      if(!exists("ke2g")) initialize()
      fileList<-paste("map",substring(pathwayList,6),".xml",sep="")     
      uGraph<-getAllPathGraph(fileList,path,"undirected",verbose)
      names(uGraph)<-pathwayList
      dGraph<-getAllPathGraph(fileList,path,"directed",verbose)
      names(uGraph)<-pathwayList
      assign("uGraph",uGraph,envir=ke2g)
      assign("dGraph",dGraph,envir=ke2g)
      assign("keggpathid2name",getPathName(),envir=ke2g)
}
#updateGraphs(mpidList)
###########################################################################
##save ke2g environment
saveKe2g<-function(file="ke2g.rda"){
      if(!exists("ke2g")) initialize()
      save(ke2g,file=file)
}
###########################################################################
##load ke2g environment
loadKe2g<-function(file="ke2g.rda"){
      if(!exists("ke2g")) initialize()
      load(ke2g,file=file)
}
###########################################################################
##load ke2g environment
getAexample<-function(k=1000){
   if(!exists("ke2g")) initialize()
   gene2path<-get("gene2path",envir=ke2g)
   if(k<=length(gene2path))
   geneList<-names(gene2path)[1:k]
   else
   geneList<-names(gene2path)
   return(geneList)
}
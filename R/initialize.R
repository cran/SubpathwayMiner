.First.lib<-function(lib, pkgname){
  #library.dynam(pkgname, pkgname, lib)
  initialize_ke2g()
  #data("hsa_ncbi-geneid");

}
########################################################################
##initialize data
initialize_ke2g<-function(){
      #assign("ke2g",new.env(parent=globalenv()),envir=.GlobalEnv)
      data("hsa_ncbi-geneid");
}

#updateOrgAndIdType("hsa","ncbi-geneid")
########################################################################
##update org and idType
updateOrgAndIdType<-function(org="hsa",idType="ncbi-geneid",
  path="ftp://ftp.genome.jp/pub/kegg/genes/organisms",verbose=TRUE){
      if(!exists("ke2g")) initialize_ke2g()
      if(verbose==TRUE){
        print("Note that the programming may be time consumming!")  
      }
	  if(verbose==TRUE){
        print("download relations between KEGG genes and current genes")
		print(paste("The current genes are genes with idType=",idType,sep=""))
	  }
	  file1<-paste(path,"/",org,"/",org,"_",idType,".list",sep="")
      #GenekeggGene<-convertFile1List(file1)
	  keggGene2gene<-read.table(file1,header=FALSE,sep = "\t", quote="\"",colClasses="character")
      assign("keggGene2gene",keggGene2gene,envir=ke2g)
      if(verbose==TRUE)
        print("download relations between KEGG genes and enzymes.")
	  file2<-paste(path,"/",org,"/",org,"_enzyme.list",sep="")
      gene2ec<-read.table(file2,header=FALSE,sep = "\t", quote="\"")
      assign("gene2ec",gene2ec,envir=ke2g)
      
      if(verbose==TRUE)
        print("download relations between KEGG genes and KOs.")
	  file3<-paste(path,"/",org,"/",org,"_ko.list",sep="")
	  gene2ko<-read.table(file3,header=FALSE,sep = "\t", quote="\"")
      assign("gene2ko",gene2ko,envir=ke2g)
	  
	  if(verbose==TRUE)
      print("download relations between KEGG gene and pathway...............")
	  file4<-paste(path,"/",org,"/",org,"_pathway.list",sep="")
	  gene2path<-read.table(file4,header=FALSE,sep = "\t", quote="\"")	  
      assign("gene2path",gene2path,envir=ke2g)
	  
	  if(verbose==TRUE)
      print("download relations between pathway identifiers and titles...............")
	  file5<-"ftp://ftp.genome.jp/pub/kegg/pathway/map_title.tab"
	  map_title<-read.table(file5,header=FALSE,sep = "\t", quote="\"",colClasses="character")	
      assign("map_title",map_title,envir=ke2g) 
	  
      assign("orgAndIdType",c(org,idType),envir=ke2g)    
}
###########################################################################
##updateGraph
updateGraphs<-function(pathwayList=getDefaultMetabolicPathway(),
         path="ftp://ftp.genome.jp/pub/kegg/release/archive/kgml/KGML_v0.6.1/map",verbose=TRUE){
      library(XML)
      if(!exists("ke2g")) initialize_ke2g()
      #assign("keggpathid2name",getPathName(),envir=ke2g)
      fileList<-paste("map",substring(pathwayList,6),".xml",sep="")     
      uGraph<-getAllPathGraph(fileList,path,"undirected",verbose)
      names(uGraph)<-pathwayList
      assign("uGraph",uGraph,envir=ke2g)   
}
#updateGraphs(mpidList)
###########################################################################
##new! updateKOGraph
updateKOGraphs<-function(pathwayList=getDefaultKOPathway(),
         path="ftp://ftp.genome.jp/pub/kegg/release/archive/kgml/KGML_v0.6.1/ko",verbose=TRUE){
      library(XML)
      if(!exists("ke2g")) initialize_ke2g()
      #assign("keggpathid2name",getPathName(),envir=ke2g)
      fileList<-paste("ko",substring(pathwayList,6),".xml",sep="")     
      KOuGraph<-getAllKOPathGraph(fileList,path,verbose)
      names(KOuGraph)<-pathwayList
      assign("KOuGraph",KOuGraph,envir=ke2g)    
}
#KOList<-list.files("e:/ko/")
#kpidList<-paste("path:",substring(KOList,3,7),sep="")
#updateKOGraphs(getDefaultKOPathway()[1])

###########################################################################
##save ke2g environment
saveKe2g<-function(file="ke2g.rda"){
      if(!exists("ke2g")) initialize_ke2g()
      save(ke2g,file=file,envir=.GlobalEnv)
}
###########################################################################
##load ke2g environment
loadKe2g<-function(file="ke2g.rda"){
      if(!exists("ke2g")) initialize_ke2g()
      load(file=file,envir=.GlobalEnv)
}
###########################################################################
##load ke2g environment
getAexample<-function(k=1000){
   if(!exists("ke2g")) initialize_ke2g()
   gene2path<-get("gene2path",envir=ke2g)
   allGene1<-getGeneFromKGene(as.character(gene2path[,1]))   
   if(k<=length(allGene1))
   geneList<-allGene1[1:k]
   else
   geneList<-allGene1

   return(geneList)
}
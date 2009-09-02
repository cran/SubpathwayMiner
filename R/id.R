#################################################################################
##get id
getkidoid<-function(org,type2,path="ftp://ftp.genome.jp/pub/kegg/genes/organisms"){
      file2<-paste(org,"_",type2,".list",sep="")
      list2<-scan(paste(path,org,file2,sep="/"),what=list(left="",right=""),sep="\t")
      result<-list()
      kgene<-sapply(strsplit(list2$left,":"),function(x) x[2])
      ogene<-sapply(strsplit(list2$right,":"),function(x) x[2])
      kgene2ogene<-ogene
      names(kgene2ogene)<-kgene
      ogene2kgene<-kgene
      names(ogene2kgene)<-ogene  
      result[[1]]<-kgene2ogene 
      result[[2]]<-ogene2kgene 
      return(result)
}
getPathName<-function(file="map_title.tab",path="ftp://ftp.genome.jp/pub/kegg/pathway_gif"){  
      list2<-scan(paste(path,file,sep="/"),what=list(left="",right=""),quote="\"",sep="\t")
      #id<-paste("path:",list2$left,sep="")  
      id<-list2$left
      name<-list2$right
      idname<-list()
      for(i in 1:length(name)){
        idname[[i]]<-name[i]
      }
      names(idname)<-id
      as.list(idname)
      return(idname)
}

#################################################################################
##merge two file and get a list of leghth=2. for example, gene2ec and ec2gene
convertFile2List<-function(file1,file2,sep="\t"){
      list1<-as.data.frame(scan(file1,what=list(left="",right=""),sep=sep))
      list2<-as.data.frame(scan(file2,what=list(left="",right=""),sep=sep))
      listm<-merge(list1,list2,by.x="left",by.y="left")
      ecList<-as.character(listm[,2])
      geneList<-sapply(strsplit(as.character(listm[,3]),":"),function(x) x[2])
      result<-list()
      uniqueList<-unique(ecList)
      result[[1]]<-sapply(uniqueList,function(y) geneList[sapply(ecList,function(x) x==y)])
      uniqueList<-unique(geneList)
      result[[2]]<-sapply(uniqueList,function(y) ecList[sapply(geneList,function(x) x==y)])
      return(result)      
}

#################################################################################
##merge two file and get a list of leghth=2. for example, gene2pathway and pathway2gene
convertFile2List1<-function(file1,file2,sep="\t"){
      list1<-as.data.frame(scan(file1,what=list(left="",right=""),sep=sep))
      list2<-as.data.frame(scan(file2,what=list(left="",right=""),sep=sep))
      listm<-merge(list1,list2,by.x="left",by.y="left")
      ecList<-paste("path",substring(as.character(listm[,2]),9),sep=":")
      geneList<-sapply(strsplit(as.character(listm[,3]),":"),function(x) x[2])
      result<-list()
      uniqueList<-unique(ecList)
      result[[1]]<-sapply(uniqueList,function(y) geneList[sapply(ecList,function(x) x==y)])
      uniqueList<-unique(geneList)
      result[[2]]<-sapply(uniqueList,function(y) ecList[sapply(geneList,function(x) x==y)])
      return(result)      
}

##############################################################################################
## merge identifiers. for example, gene2ec and ec2gene
getMerge<-function(org="hsa",type1="enzyme",type2="ncbi-geneid",path="ftp://ftp.genome.jp/pub/kegg/genes/organisms"){

      file1<-paste(org,"_",type1,".list",sep="")
      file2<-paste(org,"_",type2,".list",sep="")
      convertResult<-convertFile2List(paste(path,org,file1,sep="/"),paste(path,org,file2,sep="/"))
      return(convertResult)
}

##############################################################################################
## merge identifiers. for example, gene2pathway and pathway2gene
getMerge1<-function(org="hsa",type1="pathway",type2="ncbi-geneid",path="ftp://ftp.genome.jp/pub/kegg/genes/organisms"){

      file1<-paste(org,"_",type1,".list",sep="")
      file2<-paste(org,"_",type2,".list",sep="")
      convertResult<-convertFile2List1(paste(path,org,file1,sep="/"),paste(path,org,file2,sep="/"))
      return(convertResult)
}

##############################################################################################
## build background gene lists
buildBgGeneList<-function(org="hsa",type2="ncbi-geneid",
path="ftp://ftp.genome.jp/pub/kegg/genes/organisms"){
  file2<-paste(org,"_",type2,".list",sep="")
  list2<-as.data.frame(scan(paste(path,org,file2,sep="/"),what=list(left="",right=""),sep="\t"))
  geneList<-sapply(strsplit(as.character(list2[,2]),":"),function(x) x[2])
  uniqueList<-unique(geneList)
  return(uniqueList)
}


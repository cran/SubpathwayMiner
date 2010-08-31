printDiseAnn<-function (ann,diseCate,pathClass) 
{
if(!exists("ke2g")) initialize_ke2g();
      pathwayName<-sapply(ann,function(x) x$pathwayName)
      annGeneRatio<-sapply(ann,function(x) paste(x$annGeneNumber,x$geneNumber,sep="/"))
      annBgGeneRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
      ##<-sapply(ann,function(x) paste(x$annBgeneNumberumber,x$bgeneNumberumber,sep="/"))
      pvalue<-sapply(ann,function(x) x$pvalue)
      qvalue<-sapply(ann,function(x) x$qvalue)
      diseaseName<-sapply(ann,function(x) x$diseaseName)
      subpathwayID<-sapply(ann,function(x) x$Path)
    if(pathClass==TRUE) {pathwayClassNames<-sapply(ann,function(x) x$pathwayClassNames);}
    if(diseCate==TRUE) {diseaseCates<-sapply(ann,function(x) x$DiseaseCate);}
	  annGeneList<-sapply(ann, function(x){ paste(x$annBgGeneList,collapse=";") })
      annBgGeneList<-sapply(ann, function(x){ paste(x$annBgGeneList,collapse=";")})
      annGeneNumber<-sapply(ann,function(x) x$annGeneNumber)
	  annBgGeneNumber<-sapply(ann,function(x) x$annBgNumber)

if(diseCate==FALSE&&pathClass==TRUE) ann.data.frame<-as.data.frame(cbind(diseaseName,subpathwayID,pathwayName,annGeneRatio,annBgGeneRatio
                             ,pvalue,qvalue,annGeneNumber,annGeneList,annBgGeneNumber,annBgGeneList,pathwayClassNames));
if(pathClass==FALSE&&diseCate==TRUE) ann.data.frame<-as.data.frame(cbind(diseaseCates,diseaseName,subpathwayID,pathwayName,annGeneRatio,annBgGeneRatio
                             ,pvalue,qvalue,annGeneNumber,annGeneList,annBgGeneNumber,annBgGeneList));
if(diseCate==FALSE&&pathClass==FALSE) ann.data.frame<-as.data.frame(cbind(diseaseName,subpathwayID,pathwayName,annGeneRatio,annBgGeneRatio
                             ,pvalue,qvalue,annGeneNumber,annGeneList,annBgGeneNumber,annBgGeneList));
if(pathClass==TRUE&&diseCate==TRUE)
      ann.data.frame<-as.data.frame(cbind(diseaseCates,diseaseName,subpathwayID,pathwayName,annGeneRatio,annBgGeneRatio
                             ,pvalue,qvalue,annGeneNumber,annGeneList,annBgGeneNumber,annBgGeneList,pathwayClassNames));
    return(ann.data.frame)
}

generateNetwork<-function(inData,k=3,pvalue=0.01,geneNumber=0,diseCate=FALSE,pathClass=FALSE,exampleNumber=-1,verbose=TRUE)
{
if(!exists("ke2g")) initialize_ke2g();
pathwayClass<-get("pathwayClass",envir=ke2g);
Dise<-as.matrix(inData[1]);
Genes<-as.matrix(inData[2]);
if(ncol(inData)==2){diseCate=FALSE;}
if(diseCate==TRUE) {Cate<-as.matrix(inData[3]);}
Genes1<-matrix();
OutData<-list();
Count<-c(1);
SubGraph<-getKcSubGraph(k,graphList=getDefaultUndirectedGraph());
allDisease<-unique(Dise);
for (i in 1:length(allDisease))
{
if(i==(exampleNumber+1)){break;}
if(verbose==TRUE) print(paste("deal with the disease ",i," in ",length(allDisease)," diseases.",sep=""));
GeneNum<-which(Dise==allDisease[i]);
if (length(GeneNum)>0)
{
DiseaseAnn<-getAnn(as.character(Genes[GeneNum]),background=getDefaultBackground(),graphList=SubGraph);
if(length(DiseaseAnn)>0)
{
for(j in 1:length(DiseaseAnn))
{
if(DiseaseAnn[[j]]$pvalue<pvalue&DiseaseAnn[[j]]$annGeneNumber>geneNumber)
{
b<-list(diseaseName=allDisease[i]);
d<-list(Path=names(DiseaseAnn[j]));
e<-list(pathwayClassNames=as.character(pathwayClass[which(pathwayClass[,3]==substring(names(DiseaseAnn[j]),6,10)),1]));
if(diseCate==TRUE)
{
c<-list(DiseaseCate=as.character(Cate[GeneNum[1]]));
OutData[Count]<-list(c(c,b,d,DiseaseAnn[[j]],e));
}
else {OutData[Count]<-list(c(b,d,DiseaseAnn[[j]],e));}
Count<-Count+1;
}
}
}
}
}
pathresult<-printDiseAnn(OutData,diseCate,pathClass);
if(verbose==TRUE) print(paste("finished."));
return(pathresult)
}
#######################
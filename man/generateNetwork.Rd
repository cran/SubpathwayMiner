\name{generateNetwork}
\alias{generateNetwork}
\title{Generate disease-metabolic subpathway networks}
\description{Generate disease-metabolic subpathway networks}
\usage{
generateNetwork(inData,k=3,pvalue=0.01,geneNumber=0,
                diseCate=FALSE,pathClass=FALSE,exampleNumber=-1,verbose=TRUE)
}
\arguments{
\item{inData}{A data.frame. The disease-gene associations.}
\item{k}{A integer. The size of subpathway.}
\item{pvalue}{A numeric. The cutoff.}
\item{geneNumber}{A integer. The number of genes annotated to a subpathway.}
\item{diseCate}{A logical. It means that whether disease category will be exported.}
\item{pathClass}{A logical. It means that whether pathway class will be exported.}
\item{exampleNumber}{A integer. Number of dealing with diseases as an example of generating network. The default value is -1, dealing with all diseases.}
\item{verbose}{A logical. If TRUE, the additional diagnostics are printed}
}
\details{
The function is used to construct a disease-metabolic subpathway network (DMSPN) in which nodes represent diseases or metabolic subpathways and they are connected by an edge if genes of the disease are significantly enriched to the metabolic subpathways. We used the k-cliques subpathway identification method \code{\link{getAnn}},\code{\link{getKcSubGraph}} provided by the package to identify statistically significantly enriched disease-causing subpathways. After inputting gene sets (e.g. lung cancer) and distance parameter k, the k-clique method can mine each metabolic subpathway and then identify statistically significantly enriched subpathways. The following describes the step-by-step method for identifying metabolic subpathways. First, each metabolic pathway is converted to an undirected graph with enzymes as nodes. Enzymes in a graph are connected by an edge if their corresponding reactions have a common compound. Secondly, according to parameter k, all sub-pathways (k-cliques) in metabolic pathways can be constructed using k-cliques algorithm on each above graph. The identifier of each subpathway is given with its pathway identifier plus a sub-pathway number. For each subpathway (k-clique), distance among all enzymes within it is no greater than the parameter k (a user-defined distance). Gene sets can then be annotated to these subpathways through assigning EC numbers for them and matching them to these subpathways. Finally, the significantly enriched subpathways can be identified using hypergeometric test.To construct disease-metabolic subpathway network (DMSPN), The disease-causing gene sets need to be used to identify enriched disease-causing subpathways. Then, for each disease, the statistically significantly enriched subpathways were identified by the k-cliques subpathway identification method. Finally, the disease-metabolic subpathway network can be constructed. 
  As an example, the data used by us are obtained from the GAD(Genetic Association Database), which consist of diseases, disease-causing genes (NCBI-geneid), categories (disease classes) that the diseases belong to. The GAD is an NIH supported gene-centered public repository of human association studies examing a wide range of human diseases, including non-mendelian common diseases \url{http://geneticassociationdb.nih.gov/cgi-bin/index.cgi}.
  If \code{diseCate} is TRUE, the disease category will be exported. Note that disease category information need to be in \code{inData}.
  If \code{pathClass} is TRUE, the pathway class will be exported.

  Note that, because the function \code{generateNetwork} is time consuming, the argument \code{exampleNumber} is used to generate an example of network, which includes the given disease number \code{exampleNumber}. The default value of \code{exampleNumber}  is -1, meaning that all disease data will be treated.  

  Note that the arguments \code{pvalue} and \code{geneNumber} are important for the size of network. For example, if we set up pvalue=0.001 and geneNumber=3, then disease-subpathway associations with p>=0.001 will be ignored, and those subpathways, within which the disease gene number is less than 3, will be also ignored. 
}
\value{
A data.frame. Each row includes eleven elements at most: 'diseaseCates', 'diseaseName', 'subpathwayID', 'pathwayName', 'annGeneRatio', 'annBgGeneRatio', 'pvalue', 'qvalue', 'annGeneNumber', 'annGeneList', 'annBgGeneNumber', 'annBgGeneList', 'pathwayClassNames'. They correspond to disease category (class), disease name, subpathway id, pathway name, the ratio of the submitted genes annotated to the subpathway, the ratio of background genes annotated to the subpathway, p-value, FDR-corrected q-value, numbers of submitted genes annotated to the subpathway,
 genes (NCBI-geneid) annotated to the subpathway, numbers of background genes annotated to the subpathway, background genes annotated to the subpathway, subpathway class.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getAnn}},\code{\link{getKcSubGraph}}}
\examples{
#Use disease-gene associations from GAD.
data("inData")

#construct disease-metabolic subpathway networks. Only two diseases are treated
## as an example of using the function generateNetwork to generate networks.
DMSPN<-generateNetwork(Dise2G,exampleNumber=2)

#display 5 associations.
DMSPN[1:5,]

# save the network as .txt file.
write.table(DMSPN,"DMSPN.txt",row.names=FALSE,sep="\t")



## disease-subpathway associations with p>=0.001 are ignored, and those subpathways, 
##  within which the disease gene number is less than 3, are also ignored.
#DMSPN<-generateNetwork(Dise2G,k=4,pvalue=0.001,geneNumber=3,exampleNumber=4)

## construct networks with disease class and pathway class. 
#DMSPN<-generateNetwork(Dise2G,diseCate=TRUE,pathClass=TRUE,exampleNumber=2)


##all diseases are treated
#data("inData")
#DMSPN<-generateNetwork(Dise2G,exampleNumber=-1)

}
\keyword{file}
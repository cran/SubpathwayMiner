\name{getAnn}
\alias{getAnn}
\title{Annotate a set of genes to pathways or sub-pathways based on enzyme commission (EC)}
\description{
  Annotate a set of genes to pathways or sub-pathways based on enzyme commission (EC).
}
\usage{
  getAnn(geneList,background=getDefaultBackground(),
   order="pvalue",decreasing=FALSE,graphList=getDefaultGraph())
}
\arguments{
  \item{geneList}{A character vector of genes. }
  \item{background}{A character vector of genes used to identify the statistically significantly enriched pathways.}
  \item{order}{A character string. Should be one of "pvalue", "qvalue".}
  \item{decreasing}{A logical. Should the sort order be increasing or decreasing? }
  \item{graphList}{A list. its elements may be a \code{graph-class} or a vector of charactor.} 
}
\details{
The function can annotate a set of genes to pathways or sub-pathways and identify the statistically significantly enriched pathways.  

The default value of the argument \code{graphList} is obtained from the function \code{\link{getDefaultGraph}}. It means that genes will be annotated to all pathway. If users hope to annotate genes to metabolic pathways or sub-pathways of metabolic pathways based on enzyme commission (EC), the value of \code{graphList} should be changed with the function \code{\link{getDefaultUndirectedGraph}} or \code{\link{getKcSubGraph}}.

Before you use the function, had better use the function \code{\link{getOrgAndIdType}} to get the type of current organism and gene identifiers from  the environment variable. If the value is different from the type of organism and gene identifier in your current study, you must change them by using the function \code{\link{updateOrgAndIdType}}, \code{\link{data}} or \code{\link{loadKe2g}}.

If users don't set the argument \code{background}, the background distribution will be obtained from the whole-genome genes. Detailed information is provided in the function \code{\link{getDefaultBackground}}.
}
\value{
  A list. Each element of the list is another list. It includes eight elements: 'pathwayName', 'annGeneList', 'annGeneNumber', 'annBgNumber', 'geneNumber', 'bgNumber', 'pvalue', 'qvalue'. They correspond to pathway name, the submitted genes annotated to the pathway, numbers of submitted genes annotated to the pathway, numbers of background genes annotated to the pathway, numbers of submitted genes, numbers of background genes, p-value, and FDR-corrected q-value.


To visualize and save the results, the \code{list} can be converted to the \code{data.frame} by the function \code{\link{printAnn}}. But, note that, compared with \code{data.frame}, the \code{list} provides more information, e.g., the annotated genes are saved in the \code{list},yet not in the \code{data.frame}.

note that \code{geneList} submitted by users must be a 'character' vector.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getKOAnn}},\code{\link{getMpAnn}}, \code{\link{getKcsmpAnn}},  \code{\link{cutoffAnn}}, \code{\link{printAnn}}}

\examples{
##get a set of genes
geneList<-getAexample(k=1000)
##get the annotated results
 ann<-getAnn(geneList)
##print the annotation results to screen
 printAnn(ann)

##get the details of the results
#get the informations of the first pathway 
ann[[1]]
#get the informations of the pathway "path:00010" 
ann[["path:00010"]]
#get the genes annotated to the pathway "path:00010"
ann[["path:00010"]]$"annGeneList"

##write the annotation results to tab delimited file. 
#note that the argument col.names=NA is essential.
 geneList<-getAexample(k=1000)
 ann<-getAnn(geneList)
 result<-printAnn(ann)
 write.table(result,file="result",col.names=NA,sep="\t")

##annotate the genes to the KEGG metabolic pathways based on enzyme commission (EC)
geneList<-getAexample(k=1000)
ann<-getAnn(geneList,graphList=getDefaultUndirectedGraph())
printAnn(ann)[2:5]


##annotate the genes to the sub-pathways of metabolic 
##  pathways based on enzyme commission (EC)
geneList<-getAexample(k=1000)
subGraphList<-getKcSubGraph(k=4,graphList=getDefaultUndirectedGraph())
ann<-getAnn(geneList,graphList=subGraphList)
printAnn(ann)[2:5]


##visualize sub-pathways based on enzyme commission (EC)
plotAnn("path:00010_1",subGraphList,ann)
#go to KEGG
#plotAnn("path:00010_1",subGraphList,ann,gotoKEGG=TRUE)
}
\keyword{ file }


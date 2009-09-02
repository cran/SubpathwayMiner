\name{getKOAnn}
\alias{getKOAnn}
\title{Annotate a set of genes to sub-pathways by using KEGG Orthology (KO)}
\description{
  Annotate a set of genes to sub-pathways by using KEGG Orthology (KO).
}
\usage{
  getKOAnn(geneList,background=getDefaultBackground(),
   order="pvalue",decreasing=FALSE,graphList)
}
\arguments{
  \item{geneList}{A character vector of genes. }
  \item{background}{A character vector of genes used to identify the statistically significantly enriched pathways.}
  \item{order}{A character string. Should be one of "pvalue", "qvalue".}
  \item{decreasing}{A logical. Should the sort order be increasing or decreasing? }
  \item{graphList}{A list. its elements may be a \code{graph-class} or a vector of charactor.} 
}
\details{
The function can annotate a set of genes to sub-pathways and identify the statistically significantly enriched pathways.  

The value of \code{graphList} should be changed with the function \code{\link{getKcSubGraph}} by setting arguments \code{graphList=}\code{\link{getDefaultKOUndirectedGraph}}.

Before you use the function, had better use the function \code{\link{getOrgAndIdType}} to get the type of current organism and gene identifiers from  the environment variable. If the value is different from the type of organism and gene identifier in your current study, you must change them by using the function \code{\link{updateOrgAndIdType}}, \code{\link{data}} or \code{\link{loadKe2g}}.

If users don't set the argument \code{background}, the background distribution will be obtained from the whole-genome genes. Detailed information is provided in the function \code{\link{getDefaultBackground}}.
}
\value{
  A list. Each element of the list is another list. It includes eight elements: 'pathwayName', 'annGeneList', 'annGeneNumber', 'annBgNumber', 'geneNumber', 'bgNumber', 'pvalue', 'qvalue'. They correspond to pathway name, the submitted genes annotated to the pathway, numbers of submitted genes annotated to the pathway, numbers of background genes annotated to the pathway, numbers of submitted genes, numbers of background genes, p-value, and FDR-corrected q-value.


To visualize and save the results, the \code{list} can be converted to the \code{data.frame} by the function \code{\link{printAnn}}. But, note that, compared with \code{data.frame}, the \code{list} provides more information, e.g., the annotated genes are saved in the \code{list},yet not in the \code{data.frame}.

}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getKcSubGraph}},\code{\link{getDefaultKOUndirectedGraph}}, \code{\link{getAnn}}, \code{\link{cutoffAnn}}, \code{\link{printAnn}},\code{\link{plotKOAnn}}}

\examples{

##annotate the genes to sub-pathway by using KEGG Orthology (KO).
##get a set of genes
geneList<-getAexample(k=100)
#
subGraphList<-getKcSubGraph(k=4,graphList=getDefaultKOUndirectedGraph())
#get the annotated results
ann<-getKOAnn(geneList,graphList=subGraphList)
#print the annotation results to screen
result<-printAnn(ann)
result[2:5]

##write the annotation results to tab delimited file. 
#note that the argument col.names=NA is essential.
 write.table(result,file="result",col.names=NA,sep="\t")

}
\keyword{ file }


\name{updateGraphs}
\alias{updateGraphs}
\title{Update graphs constructed from KEGG metabolic pathways}
\description{
  Update graphs constructed from KEGG metabolic pathways.
}
\usage{
   updateGraphs(pathwayList=getDefaultMetabolicPathway(),
               path="ftp://ftp.genome.jp/pub/kegg/xml/map/",verbose=TRUE)
}
\arguments{
  \item{pathwayList}{A character vector of the KEGG metabolic pathway identifiers.}
  \item{path}{A character string.}
  \item{verbose}{A logical. If TRUE, the additional diagnostics are printed.}
}
\details{
The function is able to update the graph variable \code{\link{uGraph}} and \code{\link{dGraph}} in the environment variable of the system. 

 We construct and update the simplification version of metabolic pathways by extracting relations from the XML data of metabolic pathways. The XML data is available from the FTP site. \url{ftp://ftp.genome.jp/pub/kegg/xml/map/}


The argument \code{pathwayList} is a character vector of the KEGG metabolic pathway identifiers. The default value is provided by the return value of \code{\link{getDefaultMetabolicPathway}}. The user can change the value for constructing the subset of the graphs of the KEGG pathways.

Note that if one needs to use the updated graphs in the future, one should run the function \code{saveKe2g} to save the setting. Detailed information is provided in \code{\link{saveKe2g}} and \code{\link{loadKe2g}}.

Note that the programming is likely to be time consumming.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{uGraph}},  \code{\link{dGraph}},\code{\link{getDefaultMetabolicPathway}},\code{\link{updateOrgAndIdType}}}

\examples{
##construct the graphs list of Nucleotide Metabolism which include two pathways.
#It is considered as the subset of pathways of interest.
#updateGraphs(pathwayList=c("path:00230","path:00240"))

#get their sub-pathways
#subpathways<-getKcSubGraph(k=5,getDefaultUndirectedGraph())

#annotate genes to these sub-pathways to identify sub-pathway
##get a example of gene list
#data(gene2path)
#geneList<-names(gene2path)[1:1000]
#ann<-getAnn(geneList,graphList=subpathways)
#printAnn(ann)[2:5]

#visualize a sub-pathway of these sub-pathways
#plot(subpathways[[1]])

}
\keyword{ file }


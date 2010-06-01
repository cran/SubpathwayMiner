\name{getKcSubGraph}
\alias{getKcSubGraph}
\title{Mine sub-pathways}
\description{
  Mine sub-pathways by using the k-clique concept in social network analysis.
}
\usage{
     getKcSubGraph(k=4,graphList=getDefaultUndirectedGraph())
}
\arguments{
  \item{k}{An integer. A distance similarity parameter.}
  \item{graphList}{An undirected graph list. Detailed information is provided in \code{\link{uGraph}} or \code{\link{KOuGraph}}.}
}
\details{
 The function uses the k-clique concept in social network analysis to mine sub-pathways. In social network analysis, a k-clique in a graph is a subgraph where the distance between any two nodes is no greater than k. 

The default value of the argument \code{graphList} is obtained from the function \code{\link{getDefaultUndirectedGraph}}. It means to mine sub-pathways based on enzyme commission (EC) numbers from the version of metabolic pathways. If users hope to mine sub-pathways based on KEGG Orthology (KO), the value of \code{graphList} should be changed with the function \code{\link{getDefaultKOUndirectedGraph}}.
}
\value{
  A list of graphs.
}
\references{
Wasserman,S. and Faust,K. (1994) Social network analysis: methods and applications.  Cambridge University Press., New York, America.

Huber,W., Carey,V.J., Long,L., Falcon,S. and Gentleman,R. (2007) Graphs in molecular biology.  BMC Bioinformatics., 8, s8.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getAnn}},\code{\link{getKOAnn}},\code{\link{getKcsmpAnn}}}

\examples{
##get graph representation of metabolic pathways
graphList<-getDefaultUndirectedGraph()

#get all 4-clique subgraphs
subGraphList<-getKcSubGraph(k=4,graphList)

#display first subGraph 
library(Rgraphviz)
plot(subGraphList[[1]])



##annotate a set of genes to the sub-pathways of metabolic pathways 
##   by using enzyme commission (EC) numbers.
#geneList<-getAexample(k=1000)
#subGraphList<-getKcSubGraph(k=4)
#ann<-getAnn(geneList,graphList=subGraphList)
#printAnn(ann)

##annotate the genes to sub-pathway by using KEGG Orthology (KO).
#geneList<-getAexample(k=100)
#subGraphList<-getKcSubGraph(k=4,graphList=getDefaultKOUndirectedGraph())
#ann<-getKOAnn(geneList,graphList=subGraphList)
#printAnn(ann)[2:5]

}
\keyword{ file }


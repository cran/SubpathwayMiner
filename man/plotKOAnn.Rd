\name{plotKOAnn}
\alias{plotKOAnn}
\title{Visualize pathways or sub-pathways based on KEGG Orthology (KO)}
\description{
Visualize pathways or sub-pathways based on KEGG Orthology (KO).
}
\usage{
    plotKOAnn(pathway,graphList,ann,gotoKEGG=FALSE)
}
\arguments{
  \item{pathway}{An character string of pathway identifier, e.g., "path:00010".}
  \item{graphList}{An graph list.}
  \item{ann}{An list, e.g., the return value of the function \code{\link{getKOAnn}}.}
  \item{gotoKEGG}{An logical. Detailed informations is provided in the function \code{\link{gotoKEGG}}.}
}
\details{
The function can visualize the pathways or sub-pathways based on KEGG Orthology (KO). The red nodes in the result graph represent the enzymes which include the submitted genes.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{gotoKEGG}},\code{\link{getKOAnn}},\code{\link{printAnn}}}

\examples{
\dontrun{

##visualize sub-pathways based on KEGG Orthology (KO)
library(Rgraphviz)
geneList<-getAexample(k=1000)
subGraphList<-getKcSubGraph(k=4,graphList=getDefaultKOUndirectedGraph())
ann<-getKOAnn(geneList,graphList=subGraphList)
plotKOAnn("path:00010_1",subGraphList,ann)
#go to KEGG
plotKOAnn("path:00010_1",subGraphList,ann,gotoKEGG=TRUE)

}
}
\keyword{ file }


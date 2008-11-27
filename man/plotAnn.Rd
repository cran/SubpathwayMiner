\name{plotAnn}
\alias{plotAnn}
\title{Visualize pathways or sub-pathways}
\description{
Visualize pathways or sub-pathways.
}
\usage{
    plotAnn(pathway,graphList,ann,gotoKEGG=FALSE)
}
\arguments{
  \item{pathway}{An character string of pathway identifier, e.g., "path:00010".}
  \item{graphList}{An graph list.}
  \item{ann}{An list, e.g., the return value of the function \code{\link{getAnn}}.}
  \item{gotoKEGG}{An logical. Detailed informations is provided in the function \code{\link{gotoKEGG}}.}
}
\details{
The function can visualize the pathways or sub-pathways of metabolic pathways. The red nodes in the result graph represent the enzymes which include the submitted genes.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{gotoKEGG}},\code{\link{getAnn}},\code{\link{printAnn}}}

\examples{
##visualize metabolic pathways 
library(Rgraphviz)
geneList<-getAexample(k=1000)
ann<-getAnn(geneList,graphList=getDefaultUndirectedGraph())
plotAnn("path:00010",getDefaultUndirectedGraph(),ann)

##visualize sub-pathways of metabolic pathways
library(Rgraphviz)
geneList<-getAexample(k=1000)
subGraphList<-getKcSubGraph(k=4)
ann<-getAnn(geneList,graphList=subGraphList)
plotAnn("path:00010_1",subGraphList,ann)
#go to KEGG
plotAnn("path:00010_1",subGraphList,ann,gotoKEGG=TRUE)
}
\keyword{ file }


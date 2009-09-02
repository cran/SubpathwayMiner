\name{gotoKEGG}
\alias{gotoKEGG}
\title{Visualize pathways or sub-pathways through linking to the KEGG web site}
\description{
  Visualize pathways or sub-pathways through linking to the KEGG web site.
}
\usage{
    gotoKEGG(pathway,ann)
}
\arguments{
  \item{pathway}{A character string of pathway identifier, e.g., "path:00010".}
  \item{ann}{A list, e.g., the return value of the function \code{\link{getAnn}}.}

}
\details{
On these pictures, The red nodes represent the enzymes (or KOs) thats include the submitted genes.
.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{plotAnn}},\code{\link{getAnn}},\code{\link{printAnn}}}

\examples{
##visualize pathways through linking to the KEGG web site.
library(Rgraphviz)
geneList<-getAexample(k=100)
ann<-getAnn(geneList)
#gotoKEGG("path:00010",ann)


##visualize sub-pathways of metabolic pathways based on enzyme commission (EC) numbers through linking to the KEGG web site.
geneList<-getAexample(k=1000)
subGraphList<-getKcSubGraph(k=4)
ann<-getAnn(geneList,graphList=subGraphList)
#gotoKEGG("path:00010_1",ann)

##visualize sub-pathways sub-pathways based on KEGG Orthology (KO) identifiers through linking to the KEGG web site.
geneList<-getAexample(k=100)
subGraphList<-getKcSubGraph(k=4,graphList=getDefaultKOUndirectedGraph())
ann<-getKOAnn(geneList,graphList=subGraphList)
#gotoKEGG("path:00010_1",ann)
}
\keyword{ file }


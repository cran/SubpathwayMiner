\name{getKcsmpAnn}
\alias{getKcsmpAnn}
\title{Get the sub-pathways annotation and identification of metabolic network}
\description{
  Annotate a set of genes to sub-pathways of metabolic pathways and identify the statistically significantly enriched sub-pathways.
}
\usage{
getKcsmpAnn(geneList,background=getDefaultBackground(),k=4,
   order="pvalue",decreasing=FALSE)
}
\arguments{
  \item{geneList}{A character vector of genes. }
  \item{background}{A character vector of genes used to identify the statistically significantly enriched sub-pathways.}
  \item{k}{An integer. A distance similarity parameter.}
  \item{order}{A character string. Should be one of "pvalue", "qvalue".}
  \item{decreasing}{A logical. Should the sort order be increasing or decreasing? }
}
\details{
The function can implement the sub-pathways annotation of metabolic pathways. It is the special form of the function \code{\link{getAnn}} where the argument \code{graphList} is the return value of the function \code{\link{getKcSubGraph}}. Detailed information is provided in \code{\link{getAnn}}.

}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getMpAnn}},\code{\link{getAnn}}, \code{\link{cutoffAnn}},\code{\link{printAnn}}}

\examples{
##Annotate a set of genes to sub-pathways of metabolic pathways
# geneList<-getAexample(k=1000)
# ann<-getKcsmpAnn(geneList,k=4)
# printAnn(ann)[2:5]

}
\keyword{ file }


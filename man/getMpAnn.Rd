\name{getMpAnn}
\alias{getMpAnn}
\title{get annotation informations of metabolic pathways}
\description{
  Annotate a set of genes to metabolic pathways based on EC.
}
\usage{
  getMpAnn(geneList,background=getDefaultBackground(),
   order="pvalue",decreasing=FALSE)
}
\arguments{
  \item{geneList}{A character vector of genes.}
  \item{background}{A character vector of genes used to identify the statistically significantly enriched sub-pathways.}
  \item{order}{A character string. Should be one of "pvalue", "qvalue".}
  \item{decreasing}{A logical. Should the sort order be increasing or decreasing? }
}
\details{
The function can implement the annotation and identification of metabolic pathways. It is the special form of the function \code{\link{getAnn}} where the argument \code{graphList} is the return value of the function \code{\link{getDefaultUndirectedGraph}}. Detailed information is provided in \code{\link{getAnn}}.

}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getKcsmpAnn}},\code{\link{getAnn}}, \code{\link{cutoffAnn}},\code{\link{printAnn}}}

\examples{
\dontrun{

##annotate a set of genes to metabolic pathways
geneList<-getAexample(k=1000)
ann<-getMpAnn(geneList)
printAnn(ann)

}
}
\keyword{ file }


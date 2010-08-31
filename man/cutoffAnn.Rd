\name{cutoffAnn}
\alias{cutoffAnn}
\title{Get the statistically significantly enriched pathways}
\description{
  Get the statistically significantly enriched pathways 
  according to the results returned from \code{\link{getAnn}} or \code{\link{getKOAnn}} function.
}
\usage{
     cutoffAnn(ann,type="pvalue",operate="<=",cutoff=0.01)
}
\arguments{
  \item{ann}{ A list. The results returned from getAnn function.}
  \item{type}{A character string. Should be one of "pvalue", "qvalue".}
  \item{operate}{A character string. Should be one of "<","<=", ">", ">=".}
  \item{cutoff}{A numeric. Detailed information is provided in \code{\link{getAnn}}.}
}
\details{
 The function is used to identify the statistically signicantly enriched pathways compared with a background distribution. If users
 don't input background gene list in the function \code{\link{getAnn}} or \code{\link{getKOAnn}}, then background distribution will be obtained from the whole genome. 

}
\value{
  A list of annotation results. Detailed information is provided in \code{\link{getAnn}} or \code{\link{getKOAnn}}.
}
\references{
Strimmer,K. (2008). A unified approach to false discovery rate estimation.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getAnn}},\code{\link{getKOAnn}}, \code{\link{printAnn}}}

\examples{
\dontrun{

##get an example of gene sets
geneList<-getAexample(k=1000)

##get the annotated results
 ann<-getAnn(geneList)

##get the statistically signicantly enriched pathways according to pvalue<0.001
 cutedAnn<-cutoffAnn(ann,"pvalue","<",0.001)

##print results to screen
 printAnn(cutedAnn)

}
}
\keyword{ file }


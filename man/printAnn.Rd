\name{printAnn}
\alias{printAnn}
\title{Print the results of pathway annotation and identification}
\description{
  Print the results of pathway annotation and identification.
}
\usage{
     printAnn(ann,detail=FALSE)
}
\arguments{
  \item{ann}{ A list. The results returned from the function \code{\link{getAnn}} or \code{\link{getKOAnn}}.}
  \item{detail}{ A logical. If true, gene lists from the function \code{\link{getAnn}} or \code{\link{getKOAnn}} 
                         are converted into strings, which are used to display and write results with genes. }
}
\value{
  A data.frame of the annotation results. Its row names are pathway identifiers, 
e.g, path:00010. Columns include (pathwayName, annGeneRatio, annBgRatio, pvalue, qvalue, annGeneList, annBgGeneList). 
 
  The \code{annGeneRatio} is the ratio of the annotated genes ,e.g.,30/1000 means that 
30 genes in 1000 genes are annotated.
  The qvalue is the FDR-corrected q-value.
  annGeneList and annBgGeneList are the annotated genes and the annotated background genes.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getAnn}},\code{\link{getKOAnn}},\code{\link{cutoffAnn}}}

\examples{
##get an example of gene list
geneList<-getAexample(k=1000)

##get annotation results
 ann<-getAnn(geneList)

##print results to screen
 printAnn(ann)

##print subset of columns to visilize well
result<-printAnn(ann)
result[,2:5]

##print subset of rows to visilize well
result[1:10,]

##change print order of columns
result[,c(4,3,2,5)]


##write the annotation results to tab delimited file. 
# Notices that the argument col.names=NA is essential.
 result<-printAnn(ann)
 write.table(result,file="result",col.names=NA,sep="\t")
## write the annotation results with gene lists.
 result<-printAnn(ann,detail=TRUE)
 write.table(result,file="result",col.names=NA,sep="\t")
}
\keyword{ file }


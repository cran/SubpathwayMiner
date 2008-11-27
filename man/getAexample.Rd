\name{getAexample}
\alias{getAexample}
\title{Get a set of genes}
\description{
  Get a set of genes.
}
\usage{
 getAexample(k=1000)
}
\arguments{
  \item{k}{An integer. Set numbers of genes returned from the function.}
}

\value{
A vector of character.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\examples{
##get an example of gene sets
geneList<-getAexample(k=1000)
geneList[1:10]
}
\keyword{ file }


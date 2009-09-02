\name{getGeneFromKO}
\alias{getGeneFromKO}
\title{Get a list of genes from a list of KOs}
\description{
  Get a list of genes from a list of KOs.
}
\usage{
     getGeneFromKO(KOList)
}
\arguments{
  \item{KOList}{ A character vector of genes.}
}
\details{
Note that the result is the union of sets of genes.
}
\value{
   A character vector of genes.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getKOFromGene}}}

\examples{

## get the list of gene of two KOs. 
getGeneFromKO(c("ko:K00622","ko:K01488"))
}
\keyword{ file }


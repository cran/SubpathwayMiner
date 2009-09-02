\name{getKOFromGene}
\alias{getKOFromGene}
\title{Get a list of KOs from a list of genes}
\description{
  Get a list of KOs from a list of genes.
}
\usage{
     getKOFromGene(geneList)
}
\arguments{
  \item{geneList}{ A character vector of genes.}
}
\details{
Note that the result is the union of sets of KOs.
}
\value{
   A character vector of KOs.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getGeneFromKO}}}

\examples{
## get the list of KOs of two genes.
getKOFromGene(c("10","9","100"))

}
\keyword{ file }


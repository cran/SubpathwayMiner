\name{getPathwayFromGene}
\alias{getPathwayFromGene}
\title{Get pathwaysfrom genes}
\description{
  Get pathwaysfrom genes.
}
\usage{
     getPathwayFromGene(geneList)
}
\arguments{
  \item{geneList}{ A character vector of genes.}
}
\details{
Note that the result is the union of sets of pathways.
}
\value{
   A character vector of pathways.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getGeneFromPathway}}}

\examples{
## get the list of pathways of two genes.
getPathwayFromGene(c("5224","8802"))

}
\keyword{ file }


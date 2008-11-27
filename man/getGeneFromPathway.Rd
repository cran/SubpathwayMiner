\name{getGeneFromPathway}
\alias{getGeneFromPathway}
\title{Get genes from pathways}
\description{
  Get genes from pathways.
}
\usage{
     getGeneFromPathway(pathwayList)
}
\arguments{
  \item{pathwayList}{ A character vector of pathways.}
}
\details{
Note that the result is the union of sets of genes.
}
\value{
   A character vector of genes.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getEnzymeFromGene}}}

\examples{

## get the list of genes of two pathways. 
getGeneFromPathway(c("path:00010","path:00020"))
}
\keyword{ file }


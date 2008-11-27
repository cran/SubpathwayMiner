\name{getGeneFromEnzyme}
\alias{getGeneFromEnzyme}
\title{Get a list of genes from a list of enzymes}
\description{
  Get a list of genes from a list of enzymes.
}
\usage{
     getGeneFromEnzyme(enzymeList)
}
\arguments{
  \item{enzymeList}{ A character vector of genes.}
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

## get the list of gene of two enzymes. 
getGeneFromEnzyme(c("ec:2.7.2.3","ec:5.4.2.4"))
}
\keyword{ file }


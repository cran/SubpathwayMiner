\name{getEnzymeFromGene}
\alias{getEnzymeFromGene}
\title{Get a list of enzymes from a list of genes}
\description{
  Get a list of enzymes from a list of genes.
}
\usage{
     getEnzymeFromGene(geneList)
}
\arguments{
  \item{geneList}{ A character vector of genes.}
}
\details{
Note that the result is the union of sets of enzymes.
}
\value{
   A character vector of enzymes.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getGeneFromEnzyme}}}

\examples{
## get the list of enzymes of two genes.
getEnzymeFromGene(c("5232","5224"))

}
\keyword{ file }


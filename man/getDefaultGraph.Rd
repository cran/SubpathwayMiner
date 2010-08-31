\name{getDefaultGraph}
\alias{getDefaultGraph}
\title{Get relationship between KEGG pathways and genes}
\description{
Get relationship between KEGG pathways and genes
}
\usage{
     getDefaultGraph()
}
\details{
Get relationship between KEGG pathways and genes. The relationship is obtained from the environment variable.

Note that the "graph" is a charactor vector rather than the \code{graph-class}.
}
\value{
   A vector of charactor.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}

\examples{

pathwayId<-getDefaultGraph()

}
\keyword{ file }


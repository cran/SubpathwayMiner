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
Get relationship between KEGG pathways and genes. The relationship is obtained from the variable \code{\link{path2gene}} in the environment variable.

Note that the "graph" is a list of charactor rather than the \code{graph-class}.
}
\value{
   A list of charactor.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}

\examples{

getDefaultGraph()

}
\keyword{ file }


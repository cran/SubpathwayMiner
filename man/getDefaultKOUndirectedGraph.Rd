\name{getDefaultKOUndirectedGraph}
\alias{getDefaultKOUndirectedGraph}
\title{Get the default undirected graph based on KEGG Orthology }
\description{
Get the default undirected graph based on KEGG Orthology.
}
\usage{
     getDefaultKOUndirectedGraph
}
\details{
The default undirected graph based on KEGG Orthology is obtained from the variable \code{\link{KOuGraph}} in the environment variable.
}
\value{
   A list of \code{graph-class}.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}

\examples{

KOuGraph<-getDefaultKOUndirectedGraph()

}
\keyword{ file }


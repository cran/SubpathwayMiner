\name{getDefaultKOPathway}
\alias{getDefaultKOPathway}
\title{Get an identifier list of the default pathways based on KEGG Orthology (KO)}
\description{
  Get an identifier list of the default pathways based on KEGG Orthology (KO).
}
\usage{
     getDefaultKOPathway()
}
\details{
The identifier list of pathways based on KEGG Orthology (KO) is obtained from the variable \code{\link{kpidList}} in the environment variable.
}
\value{
   A character vector of pathway identifiers.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}

\examples{

getDefaultKOPathway()

}
\keyword{ file }


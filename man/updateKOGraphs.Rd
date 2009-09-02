\name{updateKOGraphs}
\alias{updateKOGraphs}
\title{Update graphs based on  based on KEGG Orthology (KO) constructed from KEGG pathways}
\description{
  Update graphs based on  based on KEGG Orthology (KO) constructed from KEGG metabolic pathways.
}
\usage{
updateKOGraphs(pathwayList=getDefaultKOPathway(),
         path="ftp://ftp.genome.jp/pub/kegg/xml/ko/",verbose=TRUE)
}
\arguments{
  \item{pathwayList}{A character vector of the KEGG pathway identifiers.}
  \item{path}{A character string.}
  \item{verbose}{A logical. If TRUE, the additional diagnostics are printed.}
}
\details{
The function is able to update the graph variable \code{\link{KOuGraph}} in the environment variable of the system. 

 We construct and update the simplification version by extracting relations from the XML data of pathways. The XML data is available from the FTP site. \url{ftp://ftp.genome.jp/pub/kegg/xml/ko/}

The argument \code{pathwayList} is a character vector of the KEGG metabolic pathway identifiers. The default value is provided by the return value of \code{\link{getDefaultKOPathway}}. The user can change the value for constructing the subset of the graphs of the KEGG pathways.

Note that if one needs to use the updated graphs in the future, one should run the function \code{saveKe2g} to save the setting. Detailed information is provided in \code{\link{saveKe2g}} and \code{\link{loadKe2g}}.

Note that the programming is likely to be time consumming.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{KOuGraph}},\code{\link{updateOrgAndIdType}}}

\examples{
##construct the graphs list of Nucleotide Metabolism which include two pathways.
#It is considered as the subset of pathways of interest.
#updateKOGraphs(pathwayList=c("path:00230","path:00240"))

}
\keyword{ file }


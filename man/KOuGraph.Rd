\name{KOugraph}
\docType{data}
\alias{KOuGraph}
\title{The undirected graph list based on KEGG Orthology (KO) identifiers constructed from KEGG pathways}
\description{
  Undirected graph list based on KEGG Orthology (KO) identifiers constructed from KEGG pathways.
}
\details{
For metabolic pathways, two KO identifiers are connected by an edge if there is a common compound in the KO identifiers corresponding reactions. For regulatory pathways, two KOs are connected by an edge if there are relationships between the two KOs, which can get from relation element of the XML file. The relation element specifies relationship between two KOs, which is indicated by an arrow or a line connecting two nodes in the KEGG pathways. Through the above methods to take out the relationship, all pathways are converted to an undirected graph with enzymes as nodes, which is considered as the corresponding simplification version of KEGG pathway.
}
\references{
Ogata,H., Fujibuchi,W., Goto,S. and Kanehisa,M. (2000) A heuristic graph comparison algorithm and its application to detect functionally related enzyme clusters.  
Nucleic Acids Res., 20, 4021-4028.
}
\format{A list of graph }
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{uGraph}}}
\keyword{ graphs }
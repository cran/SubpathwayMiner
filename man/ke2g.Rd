\name{ke2g}
\docType{data}
\alias{ke2g}
\alias{uGraph}
\alias{KOuGraph}
\alias{gene2path}
\alias{gene2ec}
\alias{gene2ko}
\alias{mpidList}
\alias{kpidList}
\alias{keggpathid2name}
\alias{map_title}
\title{The environment variable of the system}
\description{
  The environment variable of the system.
}
\details{
 The environment variable \code{ke2g} includes the variable \code{uGraph},\code{KOuGraph}, \code{gene2path},  \code{\link{gene2ec}}, \code{gene2ko}, \code{kpidList}, \code{mpidList}, \code{keggpathid2name}, and \code{map_title}.

Generally, A metabolic pathway can be considered as a graph with chemical compounds as nodes and enzymes as edges. We simplify metabolic pathways. Each metabolic pathway is converted to an undirected graph with enzymes as nodes. Two enzymes are connected by an edge if their corresponding reactions have a common compound. Chemical compounds are then omitted from graphs. The variable \code{uGraph} store the undirected graph list based on enzyme commission (EC) numbers constructed from KEGG metabolic pathways. Note that each graph in \code{uGraph} is constructed from metabolic pathways and one pathway can only build not more than one graph.

The variable \code{KOuGraph} store the undirected graph list based on KEGG Orthology (KO) identifiers constructed from KEGG metabolic pathways. For metabolic pathways, two KO identifiers are connected by an edge if there is a common compound in the KO identifiers corresponding reactions. For regulatory pathways, two KOs are connected by an edge if there are relationships between the two KOs, which can get from relation element of the XML file. The relation element specifies relationship between two KOs, which is indicated by an arrow or a line connecting two nodes in the KEGG pathways. Through the above methods to take out the relationship, all pathways are converted to an undirected graph with enzymes as nodes, which is considered as the corresponding simplification version of KEGG pathway.

The pathway identifier list is stored in the variable \code{kpidList}. The metabolic pathway identifier list is stored in \code{mpidList}. Detailed information is provided in \url{http://www.genome.jp/kegg/pathway.html}. The variable \code{keggpathid2name} and \code{map_title} store the data that maps KEGG pathway identifiers to names.

}
\format{An environment variable}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{loadKe2g}},\code{\link{saveKe2g}}}
\keyword{ file }
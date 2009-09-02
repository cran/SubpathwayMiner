\name{ugraph}
\docType{data}
\alias{uGraph}
\title{The undirected graph list based on enzyme commission (EC) numbers constructed from KEGG metabolic pathways}
\description{
  Undirected graph list based on enzyme commission (EC) numbers constructed from KEGG pathways.
}
\details{
 Generally, A metabolic pathway can be considered as a graph with chemical compounds as 
nodes and enzymes as edges.

We simplify metabolic pathways. Each metabolic pathway is converted to an undirected graph with enzymes as nodes. Two enzymes are connected by an edge if their corresponding reactions have a common compound. Chemical compounds are then omitted from graphs.
. 

If we consider the direction of reaction. The pathway will be a directed graph.

  Note that each graph in \code{uGraph} or \code{dGraph} is constructed from metabolic pathways
 and one pathway can only build not more than one graph.
}
\references{
Ogata,H., Fujibuchi,W., Goto,S. and Kanehisa,M. (2000) A heuristic graph comparison algorithm and its application to detect functionally related enzyme clusters.  
Nucleic Acids Res., 20, 4021-4028.
}
\format{A list of graph }
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{dGraph}},\code{\link{KOuGraph}}}
\keyword{ graphs }
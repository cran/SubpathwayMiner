\name{SubpathwayMiner-internal}
\alias{getAdjECID}
\alias{getNodes}
\alias{getPathGraph}
\alias{getAllPathGraph}
\alias{initialize_ke2g}
\alias{isKOMPath}
\alias{getKOPathGraph}
\alias{getKONodes}
\alias{getAllKOPathGraph}
\alias{getKOMPathGraph}
\alias{getAdjKOMID}
\alias{printDiseAnn}
\title{SubpathwayMiner internal functions}
\description{
  These are functions internally called by other functions
  in the package \code{SubpathwayMiner} and not meant to be called by the user.
}

\usage{
getAdjECID(ECID,top,direction)
getNodes(top)
getPathGraph(path,direction)
getAllPathGraph(fileList,path,direction,verbose=TRUE)
initialize_ke2g()
isKOMPath(top)
getKOPathGraph(top)
getKONodes(top)
getAllKOPathGraph(fileList,path,verbose=TRUE)
getKOMPathGraph(top)
getAdjKOMID(KOID,top)
printDiseAnn(ann,diseCate,pathClass)
}

\author{Chunquan Li <lcqbio@yahoo.com.cn>}

\keyword{internal}

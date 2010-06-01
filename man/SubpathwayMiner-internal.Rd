\name{SubpathwayMiner-internal}
\alias{getAdjECID}
\alias{getNodes}
\alias{getPathGraph}
\alias{getAllPathGraph}
\alias{convertFile2List}
\alias{convertFile2List1}
\alias{getMerge}
\alias{getMerge1}
\alias{buildBgGeneList}
\alias{initialize_ke2g}
\alias{getKidFromOid}
\alias{getOidFromKid}
\alias{getPathName}
\alias{getkidoid}
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
convertFile2List(file1,file2,sep="\t")
convertFile2List1(file1,file2,sep="\t")
getMerge(org="hsa",type1="enzyme",type2="ncbi-geneid",
    path="ftp://ftp.genome.jp/pub/kegg/genes/organisms")
getMerge1(org="hsa",type1="pathway",type2="ncbi-geneid",
    path="ftp://ftp.genome.jp/pub/kegg/genes/organisms")
buildBgGeneList(org="hsa",type2="ncbi-geneid",
    path="ftp://ftp.genome.jp/pub/kegg/genes/organisms")
initialize_ke2g()
getKidFromOid(oid)
getOidFromKid(kid)
getkidoid(org,type2,path="ftp://ftp.genome.jp/pub/kegg/genes/organisms")
getPathName(file="map_title.tab",
    path="ftp://ftp.genome.jp/pub/kegg/pathway_gif")
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

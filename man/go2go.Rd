\name{go2go}
\alias{getKGeneFromGene}
\alias{getGeneFromKGene}
\alias{getEnzymeFromKGene}
\alias{getKGeneFromEnzyme}
\alias{getKOFromKGene}
\alias{getKGeneFromKO}
\alias{getPathwayFromKGene}
\alias{getKGeneFromPathway}
\alias{getEnzymeFromGene}
\alias{getGeneFromEnzyme}
\alias{getKOFromGene}
\alias{getGeneFromKO}
\alias{getPathwayFromGene}
\alias{getGeneFromPathway}
\alias{getPNameFromPId}
\alias{getPIdFromPName}
\title{Map between two types of identifiers}
\description{
  Map between two types of identifiers. For example, gene identifers are converted to enzyme identifers.
}
\usage{
getKGeneFromGene(geneList)
getGeneFromKGene(keggGeneList)
getEnzymeFromKGene(keggGeneList)
getKGeneFromEnzyme(enzymeList)
getKOFromKGene(keggGeneList)
getKGeneFromKO(KOList)
getPathwayFromKGene(keggGeneList)
getKGeneFromPathway(pathwayList)
getEnzymeFromGene(geneList)
getGeneFromEnzyme(enzymeList)
getKOFromGene(geneList)
getGeneFromKO(KOList)
getPathwayFromGene(geneList)
getGeneFromPathway(pathwayList)
getPNameFromPId(PIdList)
getPIdFromPName(PNameList)
}
\arguments{
  \item{geneList}{A character vector of genes that are the current gene identifiers.}
  \item{keggGeneList}{A character vector of KEGG genes.}
  \item{enzymeList}{A character vector of enzymes.}
  \item{KOList}{A character vector of KOs.}
  \item{pathwayList}{ A character vector of pathways.}
  \item{PIdList}{A character vector of pathway identifiers.}
  \item{PNameList}{ A character vector of pathway names.}
}
\details{
Note that the argument \code{geneList} should be the current gene identifiers. The function \code{\link{getOrgAndIdType}} can get the current the type of gene identifier. The default is gene-geneid. The argument \code{keggGeneList} should be KEGG gene identifiers.

Note that the result is the union of sets of identifiers.

Note that input data should be character type. If not, the system will use the function as.character to convert data type.
}
\value{
   A character vector of identifiers.
}

\examples{
## get the list of KEGG gene identifiers of three gene(ncbi-geneid) identifiers.
getKGeneFromGene(c("1","5232","5224"))

## get the list of gene identifiers of two KEGG gene identifiers.
getGeneFromKGene(c("hsa:5232","hsa:5224"))

## get the list of enzymes of two KEGG genes.
getEnzymeFromKGene(c("hsa:5232","hsa:5224"))

## get the list of KEGG genes of two enzymes. 
getKGeneFromEnzyme(c("ec:2.7.2.3","ec:5.4.2.4"))

## get the list of KOs of two KEGG genes.
getKOFromKGene(c("hsa:5232","hsa:5224"))

## get the list of KEGG genes of KOs. 
getKGeneFromKO(c("ko:K00927","ko:K01834"))

## get the list of enzymes of two genes.
getEnzymeFromGene(c("5232","5224"))

## get the list of gene of two enzymes. 
getGeneFromEnzyme(c("ec:5.4.2.1","ec:2.7.2.3"))

## get the list of KOs of two genes.
getKOFromGene(c("5232","5224"))

## get the list of gene of two KOs.
getGeneFromKO(c("ko:K01834","ko:K00927"))

## get the list of pathways of two genes.
getPathwayFromGene(c("5224","8802"))

## get the list of genes of two pathways. 
getGeneFromPathway(c("path:00010","path:00020"))

## get pathway name. 
getPNameFromPId(c("path:00010","path:00020"))

}
\keyword{ file }


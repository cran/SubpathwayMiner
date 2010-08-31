\name{updateOrgAndIdType}
\alias{updateOrgAndIdType}
\title{Update the organism and the type of gene identifier}
\description{
  Update the organism and the type of gene identifiers.
}
\usage{
     updateOrgAndIdType(org="hsa",idType="ncbi-geneid",
   path="ftp://ftp.genome.jp/pub/kegg/genes/organisms",verbose=TRUE)
}
\arguments{
  \item{org}{A character string. The abbreviation of a genome name.}
  \item{idType}{A character string. The type of gene identifier.}
  \item{path}{A character string.}
  \item{verbose}{A logical. If TRUE, the additional diagnostics are printed.}
}
\details{
The existing tools mainly use DBMS (data base management system) to store all data relative to analysis of pathways and the update process of the data is transparent to users, which means that the annotation results users get from these tools may become outdated. We don't use DBMS to store data. We present a new method that enables users to update data by themselves. Users are firstly required to set organism and type of gene identifier before annotateing genes to the pathways. According to the setting, the system can download all data relative to analysis of pathways in the certain organism, and then treat and store them in an environment variable in R. Through the method the system can synchronize the data with the KEGG databases and support almost all organisms and cross reference identifiers in the KEGG GENE database. 

The function is able to update the variables in the environment variable \code{\link{ke2g}}, including \code{gene2ec}, \code{gene2ko}, \code{gene2path}, etc. 

Note that if the user don't run the function \code{updateOrgAndIdType}, \code{\link{loadKe2g}} after starting up R-system and loading the package of the system, then the default value of the argument \code{org} and \code{idType} is "hsa" (human) and "ncbi-geneid" (Entrez gene identifiers). The user can get the information from the return value of the funciton \code{\link{getOrgAndIdType}}.

The argument \code{org} must be the abbreviation of a genome name. For example, the ath is the abbreviation of Arabidopsis thaliana. The osa is the abbreviation of Oryza sativa japonica. The dme is the abbreviation of Drosophila melanogaster. The eco is the abbreviation of Escherichia coli K-12 MG1655. The hsa is the abbreviation of Homo sapiens. The sce is the abbreviation of Saccharomyces cerevisiae. The mmu is the abbreviation of Mus musculus. Cel is the abbreviation of Caenorhabditis elegans. Detailed information is provided in \url{http://www.genome.jp/kegg/catalog/org_list.html}.

In old versions (<3.1), the environment variable of some organisms were provided. In the new version, These data are provided not longer. One reason is that it is not very necessary to store outdated data since we can use new method to store and update data, which significantly reduce time in obtaining up-to-date data. Another is that saving size of the software is beneficial.

The argument \code{idType} is a character string of the type of identifier. The system supports most KEGG cross-reference identifiers such as Entrez gene IDs (idType="ncbi-geneid"), NCBI gi numbers (idType="ncbi-gi"), UniProt accession numbers (idType="uniprot"), etc. Detailed information is provided in \url{ftp://ftp.genome.jp/pub/kegg/genes/organisms}. For example, because a file name in "hsa" file directory is "hsa ensembl-hsa.list", idType="ensembl-hsa" is available as the input identifier type. Note that the \code{idType} is relative to the genome. Different genomes may support different \code{idType}. For example, "sgd-sce" is supported by yeast genome. however, it is not supported by human genome.

The argument \code{path} is the path of file directory of the organism cross-reference identifiers. The default value is "ftp://ftp.genome.jp/pub/kegg/genes/organisms". The setting ensure that the user is able to obain the updated data from the KEGG FTP site. Of course, the user can also download the organisms data of interest from FTP site and change path to the data file for implementing the local update.

Note that the programming is time consumming. For Sovling the problem, see \code{\link{saveKe2g}} and \code{\link{loadKe2g}}.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{getAnn}}, \code{\link{updateGraphs}},\code{\link{getOrgAndIdType}}}

\examples{
\dontrun{
##update organism and the type of gene identifiers
getOrgAndIdType()
updateOrgAndIdType("sce","sgd-sce")
getOrgAndIdType()

}
}
\keyword{ file }


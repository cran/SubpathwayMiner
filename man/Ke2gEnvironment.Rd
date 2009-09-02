\name{Ke2gEnvironment}
\alias{saveKe2g}
\alias{loadKe2g}
\title{ke2g environment function}
\description{
  Load and save the  environment variable ke2g of the system.
}
\usage{
    saveKe2g(file="ke2g.rda")
    loadKe2g(file="ke2g.rda")
}
\arguments{
  \item{file}{ A character string.}
}
\details{
  The functions are used to load or save the environment variable \code{\link{ke2g}} of the system. If one has changed the environment variable with some functions, e.g., \code{\link{updateOrgAndIdType}} or \code{\link{updateGraphs}} or \code{\link{updateKOGraphs}} and hope to use the setting in the future, then the functions can solve the problem. For example, one uses the function \code{saveKe2g} to save the environment variable. When one needs to use the setting next time, one can use the function \code{\link{loadKe2g}} to load the last environment variable.

  The functions implement the localization of the system. Therefore, the system not only implements  
the most up-to-date annotations but also the localized annotations. The functions are very important for the user to frequently annotate genes 
with the different genomes and the type of identifiers. It can increase largely the running speed because one can update 
the data one time only, and then repeatedly use it.
 
Note that if and only if the function \code{loadKe2g} runs after loading the \code{SubpathwayMiner} package, it will be in effect.



In addition, the environment variables of organisms with well annotated genomes are provided by the system and users can use the function \code{\link{data}} to load them.

}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{updateOrgAndIdType}},\code{\link{updateGraphs}},\code{\link{updateKOGraphs}}}

\examples{
##for example 
#save two ke2g environment about yeast. 
#Note that the data is saved to the working directory 
#and its name is sce_ncbi-geneid.rda.
# updateOrgAndIdType("sce","ncbi-geneid")
# saveKe2g("sce_ncbi-geneid.rda")

#shut down the R system
#start up the R system
#library(SubpathwayMiner)

#when the user annotate yeast genes. Note that you need to 
#working directory to the directory of the data file.
#loadKe2g("sce_ncbi-geneid.rda") 

#The code below is used to load the environment variables provided by the system. 
#the type of gene identifier is ncbi-geneid.
#data("sce_ncbi-geneid")

}
\keyword{ file }


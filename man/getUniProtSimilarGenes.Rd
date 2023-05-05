\name{getUniProtSimilarGenes}
\alias{getUniProtSimilarGenes}
\title{Get similar UniProt genes by a given identity percentage}
\description{
Retrieves the cluster of genes with a certain percentage of identity where a given UniProt gene belongs.}
\usage{
getUniProtSimilarGenes(upId, clusterIdentity = '1.0', clusterNames = FALSE)
}
\arguments{
    \item{upId}{valid UniProt ID}
    \item{clusterIdentity}{identity of the cluster to be retrieved. Default value '1.0' returns
    the 100\% identical genes cluster, values '0.9' and '0.5', the 90\% and 50\% identical genes clusters,
    respectively}
    \item{clusterNames}{logical value (TRUE or FALSE) indicating if cluster names are retrieved}
}
\value{
Depending on parameter 'clusterNames': if TRUE it will return a named list with cluster names as keys and
a vector of gene IDs as values, empty list if no clusters were found; if FALSE, returns a vector of UniProt
IDs of the cluster genes, empty vector if no clusters or UniProt IDs were found.
}
\examples{
# Get UniProt IDs of 50\% identical genes to UniProt gene 'G0L217'
if(interactive()){
    getUniProtSimilarGenes('G0L217', clusterIdentity = '0.5')
}
}
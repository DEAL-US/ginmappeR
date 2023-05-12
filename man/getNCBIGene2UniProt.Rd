\name{getNCBIGene2UniProt}
\alias{getNCBIGene2UniProt}
\title{Translates a NCBI Gene database ID to UniProt database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the UniProt database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBIGene2UniProt(ncbiId, exhaustiveMapping = FALSE,
byIdenticalProteins = TRUE)
}
\arguments{
  \item{ncbiId}{valid NCBI Gene database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
  \item{byIdenticalProteins}{logical value (TRUE or FALSE) indicating to
  retrieve translations of NCBI Identical Proteins to the given identifier}
}
\value{
A vector of UniProt IDs, empty if no translations were found.
}
\examples{
# Translate NCBI Gene ID '76524190' to UniProt
if(interactive()){
    getNCBIGene2UniProt('76524190')
}
}

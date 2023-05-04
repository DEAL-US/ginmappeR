\name{getNCBIProtein2UniProt}
\alias{getNCBIProtein2UniProt}
\title{Translates an NCBI Protein database ID to UniProt database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database ID to
the UniProt database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBIProtein2UniProt(ncbiId, exhaustiveMapping = FALSE, byIdenticalProteins = TRUE)
}
\arguments{
  \item{ncbiId}{valid NCBI Protein database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve all possible translations
  or only the first and fastest one}
  \item{byIdenticalProteins}{logical value (TRUE or FALSE) indicating to retrieve translations of
  NCBI Identical Proteins to the given identifier.}
}
\value{
A vector of UniProt IDs, empty if no translations were found.
}
\examples{
# Translate NCBI Protein ID 'WP_010896559.1' to UniProt
if(interactive()){
    getNCBIProtein2UniProt('WP_010896559.1')
}
}

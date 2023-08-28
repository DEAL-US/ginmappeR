\name{getNCBIProtein2UniProt}
\alias{getNCBIProtein2UniProt}
\title{Translates a NCBI Protein database ID to UniProt database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database
ID to the UniProt database.}
\usage{
getNCBIProtein2UniProt(ncbiId, exhaustiveMapping = FALSE,
detailedMapping = FALSE, byIdenticalProteins = TRUE)
}
\arguments{
  \item{ncbiId}{valid NCBI Protein database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
  \item{detailedMapping}{logical value (TRUE or FALSE) indicating whether to
  return a vector of UniProt IDs if FALSE, or a named list if TRUE. This named
  list details the origin of the translated IDs, dividing them in 'DT', if it
  is a Direct Translation, or '1.0' (100\% identity), if it has been retrieved
  through Identical Proteins database. See examples below for further 
  understanding}
  \item{byIdenticalProteins}{logical value (TRUE or FALSE) indicating to
  retrieve translations of NCBI Identical Proteins to the given identifier}
}
\value{
Depending on parameter 'detailedMapping': if TRUE it will return a named list
dividing UniProt IDs in 'DT', if it is a Direct Translation, or '1.0', if it
has been retrieved through NCBI Identical Protein database , or empty list if
no translations were found; if FALSE, returns a vector of UniProt IDs, empty
vector if no translations were found.}
\examples{
# Translate NCBI Protein ID 'WP_010896559.1' to UniProt
if(interactive()){
    getNCBIProtein2UniProt('WP_010896559.1')
}
# Translate NCBI Protein ID 'WP_010896559.1' to UniProt obtaining all possible
# translations
if(interactive()){
    getNCBIProtein2UniProt('WP_010896559.1',exhaustiveMapping = TRUE,
    detailedMapping = TRUE, byIdenticalProteins = FALSE)}
}

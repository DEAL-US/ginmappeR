\name{getNCBIProtein2CARD}
\alias{getNCBIProtein2CARD}
\title{Translates a NCBI Protein database ID to CARD database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database
ID to the CARD database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBIProtein2CARD(ncbiId, exhaustiveMapping = FALSE)
}
\arguments{
  \item{ncbiId}{valid NCBI Protein database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
}
\value{
A vector of CARD IDs, empty if no translations were found.
If an error occurs, returns NULL.
}
\examples{
# Translate NCBI Protein ID 'AAK64581' to CARD
getNCBIProtein2CARD('AAK64581')

# Translate NCBI Protein IDs 'AAK64581' and 'WP_010896559.1' to CARD
getNCBIProtein2CARD(c('AAK64581', 'WP_010896559.1'))
}

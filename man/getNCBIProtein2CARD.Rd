\name{getNCBIProtein2CARD}
\alias{getNCBIProtein2CARD}
\title{Translates a NCBI Protein database ID to CARD database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database ID to
the CARD database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBIProtein2CARD(ncbiId, exhaustiveMapping = FALSE)
}
\arguments{
  \item{ncbiId}{valid NCBI Protein database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve all possible translations
  or only the first and fastest one}
}
\value{
A vector of CARD IDs, empty if no translations were found.
}
\examples{
# Translate NCBI Protein ID 'AAK64581' to CARD
if(interactive()){
    getNCBIProtein2CARD('AAK64581')
}
# Translate NCBI Protein ID 'CAC81324.1' to CARD
if(interactive()){
    getNCBIProtein2CARD('CAC81324.1')
    }
}

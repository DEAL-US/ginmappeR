\name{getNCBIProtein2CARD}
\alias{getNCBIProtein2CARD}
\title{Translates a NCBI Protein database ID to CARD database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database
ID to the CARD database.}
\usage{
getNCBIProtein2CARD(ncbiId, exhaustiveMapping = FALSE)
}
\arguments{
  \item{ncbiId}{valid NCBI Protein database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (\code{TRUE} or \code{FALSE})
  indicating to retrieve all possible translations or only the first and
  fastest one.}
}
\value{
A character vector of the same length as \code{ncbiId}. If parameter
\code{exhaustiveMapping} is \code{TRUE}, it will return a list instead of a
character vector to avoid mixing the identifiers.
}
\examples{
# Translate NCBI Protein ID 'AAK64581' to CARD
getNCBIProtein2CARD('AAK64581')
\donttest{
# Translate NCBI Protein IDs 'AAK64581' and 'WP_010896559.1' to CARD
getNCBIProtein2CARD(c('AAK64581', 'WP_010896559.1'))

# Get all possible translations of IDs 'AAK64581', 'WP_010896559.1'
getNCBIProtein2CARD(c('AAK64581', 'WP_010896559.1'), exhaustiveMapping = TRUE)
}
}

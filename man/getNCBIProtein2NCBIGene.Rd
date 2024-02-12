\name{getNCBIProtein2NCBIGene}
\alias{getNCBIProtein2NCBIGene}
\title{Translates a NCBI Protein database ID to NCBI Gene database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database
ID to the NCBI Gene database.}
\usage{
getNCBIProtein2NCBIGene(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Protein database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (\code{TRUE} or \code{FALSE})
  indicating to retrieve all possible translations or only the first and
  fastest one.}
}
\value{
A character vector of the same length as \code{id}. If parameter
\code{exhaustiveMapping} is \code{TRUE}, it will return a list instead of a
character vector to avoid mixing the identifiers.
}
\examples{
# Translate NCBI Protein ID 'CAA79696' to NCBI Gene
getNCBIProtein2NCBIGene('CAA79696')
\donttest{
# Translate NCBI Protein IDs 'CAA79696' and 'WP_001082319' to NCBI Gene
getNCBIProtein2NCBIGene(c('CAA79696', 'WP_001082319'))

# Get all possible translations of IDs 'CAA79696', 'WP_001082319'
getNCBIProtein2NCBIGene(c('CAA79696', 'WP_001082319'), exhaustiveMapping=TRUE)
}
}

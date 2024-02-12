\name{getNCBIGene2CARD}
\alias{getNCBIGene2CARD}
\title{Translates a NCBI Gene database ID to CARD database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the CARD database.}
\usage{
getNCBIGene2CARD(ncbiId, exhaustiveMapping = FALSE)
}
\arguments{
  \item{ncbiId}{valid NCBI Gene database ID or a vector of them}
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
\donttest{
# Translate NCBI Gene ID '3510143' to CARD
getNCBIGene2CARD('3510143')

# Translate NCBI Gene IDs '3510143' and '29426913' to CARD
getNCBIGene2CARD(c('3510143', '29426913'))

# Get all possible translations of IDs '3510143' and '29426913'
getNCBIGene2CARD(c('3510143', '29426913'), exhaustiveMapping = TRUE)
}
}

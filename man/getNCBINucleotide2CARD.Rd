\name{getNCBINucleotide2CARD}
\alias{getNCBINucleotide2CARD}
\title{Translates a NCBI Nucleotide database ID to CARD database}
\description{
Retrieves the ID translation of an specific instance of NCBI Nucleotide database
ID to the CARD database.}
\usage{
getNCBINucleotide2CARD(ncbiId, exhaustiveMapping = FALSE)
}
\arguments{
  \item{ncbiId}{valid NCBI Nucleotide database ID or a vector of them}
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
# Translate NCBI Nucleotide ID 'AY034138' to CARD
getNCBINucleotide2CARD('AY034138')

# Translate NCBI Nucleotide IDs 'AY034138' and 'X54723.1' to CARD
getNCBINucleotide2CARD(c('AY034138', 'X54723.1'))

# Get all possible translations of IDs 'AY034138', 'X54723.1'
getNCBINucleotide2CARD(c('AY034138', 'X54723.1'), exhaustiveMapping = TRUE)
}
}

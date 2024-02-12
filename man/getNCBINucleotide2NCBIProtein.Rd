\name{getNCBINucleotide2NCBIProtein}
\alias{getNCBINucleotide2NCBIProtein}
\title{Translates a NCBI Nucleotide database ID to NCBI Protein database}
\description{
Retrieves the ID translation of an specific instance of NCBI Nucleotide database
ID to the NCBI Protein database.}
\usage{
getNCBINucleotide2NCBIProtein(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Nucleotide database ID or a vector of them}
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
# Translate NCBI Nucleotide ID 'AY536519' to NCBI Protein
getNCBINucleotide2NCBIProtein('AY536519')
\donttest{
# Translate NCBI Nucleotide IDs 'AY536519' and 'X54723.1' to NCBI Protein
getNCBINucleotide2NCBIProtein(c('AY536519', 'X54723.1'))

# Get all possible translations of IDs 'AY536519', 'X54723.1'
getNCBINucleotide2NCBIProtein(c('AY536519', 'X54723.1'),
exhaustiveMapping = TRUE)
}
}

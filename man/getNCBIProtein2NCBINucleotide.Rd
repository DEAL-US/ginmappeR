\name{getNCBIProtein2NCBINucleotide}
\alias{getNCBIProtein2NCBINucleotide}
\title{Translates a NCBI Protein database ID to NCBI Nucleotide database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database
ID to the NCBI Nucleotide database. Note that, due to the nature of
NCBI Nucleotide database, a complete genome is often retrieved.}
\usage{
getNCBIProtein2NCBINucleotide(id, exhaustiveMapping=FALSE)
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
\donttest{
# Translate NCBI Protein ID 'AAS48620' to NCBI Nucleotide
getNCBIProtein2NCBINucleotide('AAS48620')

# Translate NCBI Protein IDs 'AAS48620' and 'CAA38525' to NCBI Nucleotide
getNCBIProtein2NCBINucleotide(c('AAS48620', 'CAA38525'))

# Get all possible translations of IDs 'AAS48620', 'CAA38525'
getNCBIProtein2NCBINucleotide(c('AAS48620', 'CAA38525'), exhaustiveMapping=TRUE)
}
}

\name{getNCBIGene2NCBINucleotide}
\alias{getNCBIGene2NCBINucleotide}
\title{Translates a NCBI Gene database ID to NCBI Nucleotide database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the NCBI Nucleotide database. Note that, due to the nature of
NCBI Nucleotide database, a complete genome is often retrieved.}
\usage{
getNCBIGene2NCBINucleotide(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Gene database ID or a vector of them}
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
# Translate NCBI Gene ID '76524190' to NCBI Nucleotide
getNCBIGene2NCBINucleotide('76524190')
\donttest{
# Translate NCBI Gene IDs '76524190' and '1272' to NCBI Nucleotide
getNCBIGene2NCBINucleotide(c('76524190', '1272'))

# Get all possible translations of IDs '76524190' and '1272'
getNCBIGene2NCBINucleotide(c('76524190', '1272'), exhaustiveMapping = TRUE)
}
}

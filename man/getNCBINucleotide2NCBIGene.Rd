\name{getNCBINucleotide2NCBIGene}
\alias{getNCBINucleotide2NCBIGene}
\title{Translates a NCBI Nucleotide database ID to NCBI Gene database}
\description{
Retrieves the ID translation of an specific instance of NCBI Nucleotide database
ID to the NCBI Gene database.}
\usage{
getNCBINucleotide2NCBIGene(id, exhaustiveMapping=FALSE)
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
\donttest{
# Translate NCBI Nucleotide ID 'Z21488' to NCBI Gene
getNCBINucleotide2NCBIGene('Z21488')

# Get all possible translations of ID 'Z21488'
getNCBINucleotide2NCBIGene(c('Z21488'), exhaustiveMapping = TRUE)
}
}

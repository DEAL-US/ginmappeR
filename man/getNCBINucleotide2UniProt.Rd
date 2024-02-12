\name{getNCBINucleotide2UniProt}
\alias{getNCBINucleotide2UniProt}
\title{Translates a NCBI Nucleotide database ID to UniProt database}
\description{
Retrieves the ID translation of an specific instance of NCBI Nucleotide database
ID to the UniProt database.}
\usage{
getNCBINucleotide2UniProt(ncbiId, exhaustiveMapping = FALSE,
detailedMapping = FALSE, byIdenticalProteins = TRUE)
}
\arguments{
  \item{ncbiId}{valid NCBI Nucleotide database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (\code{TRUE} or \code{FALSE})
  indicating to retrieve all possible translations or only the first and
  fastest one.}
  \item{detailedMapping}{logical value (\code{TRUE} or \code{FALSE}) indicating
  to return a vector of UniProt IDs if \code{FALSE}, or a named vector if
  \code{TRUE}. This named vector details the origin of the translated IDs,
  dividing them in \code{DT}, if it is a Direct Translation, or \code{1.0},
  \code{0.9} or \code{0.5}, if it has been retrieved through any of the
  identical gene clusters.
  }
  \item{byIdenticalProteins}{logical value (\code{TRUE} or \code{FALSE})
  indicating to translate also through NCBI Identical Proteins. \code{TRUE}
  value could mean higher translation time and \code{FALSE}, fewer chances to
  translate the ID}
}
\value{
A character vector of the same length as \code{ncbiId}. If parameter
\code{exhaustiveMapping} is \code{TRUE}, it will return a list instead of a
character vector to avoid mixing the identifiers.
}
\examples{
# Translate NCBI Nucleotide ID 'AY536519' to UniProt
getNCBINucleotide2UniProt('AY536519')
\donttest{
# Translate NCBI Nucleotide IDs 'AY536519' and 'X54723.1' to UniProt
getNCBINucleotide2UniProt(c('AY536519', 'X54723.1'))

# Get all possible translations of IDs 'AY536519', 'X54723.1'
# with detailed origin
getNCBINucleotide2UniProt(c('AY536519', 'X54723.1'), detailedMapping = TRUE,
exhaustiveMapping = TRUE)
}
}

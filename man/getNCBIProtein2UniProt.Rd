\name{getNCBIProtein2UniProt}
\alias{getNCBIProtein2UniProt}
\title{Translates a NCBI Protein database ID to UniProt database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database
ID to the UniProt database.}
\usage{
getNCBIProtein2UniProt(ncbiId, exhaustiveMapping = FALSE,
detailedMapping = FALSE, byIdenticalProteins = TRUE)
}
\arguments{
  \item{ncbiId}{valid NCBI Protein database ID or a vector of them}
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
# Translate NCBI Protein ID 'WP_010896559.1' to UniProt
getNCBIProtein2UniProt('WP_010896559.1')
\donttest{
# Translate NCBI Protein IDs 'AAS48620' and 'WP_010896559' to UniProt
getNCBIProtein2UniProt(c('AAS48620', 'WP_010896559'))

# Get all possible translations of IDs 'AAS48620', 'WP_010896559'
# with detailed origin
getNCBIProtein2UniProt(c('AAS48620', 'WP_010896559'), detailedMapping = TRUE,
exhaustiveMapping = TRUE)
}
}

\name{getUniProt2CARD}
\alias{getUniProt2CARD}
\title{Translates an UniProt database ID to CARD database}
\description{
Retrieves the ID translation of an specific instance of UniProt database ID to
the CARD database.}
\usage{
getUniProt2CARD(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE,
bySimilarGenes = TRUE)
}
\arguments{
  \item{upId}{valid UniProt database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (\code{TRUE} or \code{FALSE})
  indicating to retrieve all possible translations or only the first and
  fastest one.}
  \item{detailedMapping}{logical value (\code{TRUE} or \code{FALSE}) indicating
  to return a vector of CARD IDs if \code{FALSE}, or a named vector if
  \code{TRUE}. This named vector details the origin of the translated IDs,
  dividing them in \code{DT}, if it is a Direct Translation, or \code{1.0},
  \code{0.9} or \code{0.5}, if it has been retrieved through any of the
  identical gene clusters.
  }
  \item{bySimilarGenes}{logical value (\code{TRUE} or \code{FALSE}) indicating
  to translate also through similar UniProt genes (100\%, 90\% or 50\% of
  identity). Using \code{TRUE} value may take a long time, only advised if
  there is not a direct translation possible.
  }
}
\value{
A character vector of the same length as \code{upId}. If parameter
\code{exhaustiveMapping} is \code{TRUE}, it will return a list instead of a
character vector to avoid mixing the identifiers.
}
\examples{
# Translate UniProt ID 'A0A1S7BGS4' to CARD
getUniProt2CARD('A0A1S7BGS4')
\donttest{
# Translate UniProt IDs 'A0A1S7BGS4' and 'A0A6H2TXZ6' to CARD
getUniProt2CARD(c('A0A1S7BGS4', 'A0A6H2TXZ6'))

# Get all possible translations of IDs 'A0A1S7BGS4', 'A0A6H2TXZ6'
# with detailed origin
getUniProt2CARD(c('A0A1S7BGS4', 'A0A6H2TXZ6'), detailedMapping = TRUE,
exhaustiveMapping = TRUE)
}
}

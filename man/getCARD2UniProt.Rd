\name{getCARD2UniProt}
\alias{getCARD2UniProt}
\title{Translates a CARD database ID to UniProt database}
\description{
Retrieves the ID translation of an specific instance of CARD database ID to
the UniProt database.}
\usage{
getCARD2UniProt(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE)
}
\arguments{
  \item{cardId}{valid CARD database ARO ID or a vector of them}
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
}
\value{
A character vector of the same length as \code{cardId}. If parameter
\code{exhaustiveMapping} is \code{TRUE}, it will return a list instead of a
character vector to avoid mixing the identifiers.
}
\examples{
\donttest{
# Translate CARD ID '3002867' to UniProt
getCARD2UniProt('3002867')

# Translate CARD IDs '3002867' and '3002535' to UniProt
getCARD2UniProt(c('3002867', '3002535'))

# Get all possible translations of IDs '3002867', '3002535' with detailed origin
getCARD2UniProt(c('3002867', '3002535'), detailedMapping = TRUE,
exhaustiveMapping = TRUE)
}
}

\name{getCARD2NCBIGene}
\alias{getCARD2NCBIGene}
\title{Translates a CARD database ID to NCBI Gene database}
\description{
Retrieves the ID translation of an specific instance of CARD database ID to
the NCBI Gene database.}
\usage{
getCARD2NCBIGene(cardId, exhaustiveMapping = FALSE)
}
\arguments{
  \item{cardId}{valid CARD database ARO ID or a vector of them}
  \item{exhaustiveMapping}{logical value (\code{TRUE} or \code{FALSE})
  indicating to retrieve all possible translations or only the first and
  fastest one.}
}
\value{
A character vector of the same length as \code{cardId}. If parameter
\code{exhaustiveMapping} is \code{TRUE}, it will return a list instead of a
character vector to avoid mixing the identifiers.
}
\examples{
# Translate CARD ID '3002525' to NCBI Gene
getCARD2NCBIGene('3002525')
\donttest{
# Translate CARD IDs '3002524' and '3002525' to NCBI Gene
getCARD2NCBIGene(c('3002524', '3002525'))

# Get all possible translations of IDs '3005061' and '3002524'
getCARD2NCBIGene(c('3005061', '3002524'), exhaustiveMapping = TRUE)
}
}

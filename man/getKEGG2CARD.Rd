\name{getKEGG2CARD}
\alias{getKEGG2CARD}
\title{Translates a KEGG database ID to CARD database}
\description{
Retrieves the ID translation of an specific instance of KEGG Genes database
ID to the CARD database.}
\usage{
getKEGG2CARD(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE,
bySimilarGenes = TRUE)
}
\arguments{
  \item{keggId}{valid KEGG Genes database ID or a vector of them}
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
A character vector of the same length as \code{keggId}. If parameter
\code{exhaustiveMapping} is \code{TRUE}, it will return a list instead of a
character vector to avoid mixing the identifiers.
}
\examples{
# Translate KEGG ID 'ag:ACC85616' to CARD
getKEGG2CARD('ag:ACC85616')
\donttest{
# Translate KEGG IDs 'ag:ACC85616' and 'ag:CAJ47134' to CARD
getKEGG2CARD(c('ag:ACC85616', 'ag:CAJ47134'))

# Get all possible translations of IDs 'ag:ACC85616', 'ag:CAJ47134'
# with detailed origin
getKEGG2CARD(c('ag:ACC85616', 'ag:CAJ47134'), detailedMapping = TRUE,
exhaustiveMapping = TRUE)
}
}

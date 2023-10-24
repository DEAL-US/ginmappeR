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
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
  \item{detailedMapping}{logical value (TRUE or FALSE) indicating whether to
  return a vector of UniProt IDs if FALSE, or a named list if TRUE. This named
  list details the origin of the translated IDs, dividing them in 'DT', if it
  is a Direct Translation, or '1.0', if it has been retrieved through Identical
  Proteins database. See examples below for further understanding}
}
\value{
Depending on parameter 'detailedMapping': if TRUE it will return a named list
dividing UniProt IDs in 'DT', if it is a Direct Translation (100\% of identity),
or '1.0', if it has been retrieved through NCBI Identical Protein database, or
empty list if no translations were found; if FALSE, returns a vector of UniProt
IDs, empty vector if no translations were found.}
\examples{
# Translate CARD ID '3002867' to UniProt
getCARD2UniProt('3002867')
}

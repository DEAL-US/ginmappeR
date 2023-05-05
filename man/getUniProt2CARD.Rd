\name{getUniProt2CARD}
\alias{getUniProt2CARD}
\title{Translates an UniProt database ID to CARD database}
\description{
Retrieves the ID translation of an specific instance of UniProt database ID to
the CARD database.}
\usage{
getUniProt2CARD(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE)
}
\arguments{
  \item{upId}{valid UniProt database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve all possible translations
  or only the first and fastest one}
  \item{detailedMapping}{logical value (TRUE or FALSE) indicating to return a vector of CARD IDs if FALSE,
  or a named list if TRUE. This named list details the origin of the translated IDs, dividing them in 'DT',
  if it is a Direct Translation, or '1.0', '0.9' or '0.5', if it has been retrieved through any of the identical
  gene clusters. See examples below for further understanding
  }
  \item{bySimilarGenes}{logical value (TRUE or FALSE) indicating to translate also through similar
  UniProt genes (100\%, 90\% or 50\% of identity). Using TRUE value may take a long time, only advised if
  there is not a direct translation possible
  }
}
\value{
Depending on parameter 'detailedMapping': if TRUE it will return a named list dividing CARD IDs in 'DT', if
it is a Direct Translation, or '1.0', '0.9' or '0.5', if it has been retrieved through any of the identical
gene clusters, or empty list if no translations were found; if FALSE, returns a vector of CARD IDs, empty
vector if no translations were found.
}
\examples{
# Translate UniProt ID 'A0A1S7BGS4' to CARD
if(interactive()){
    getUniProt2CARD('A0A1S7BGS4')
}
# Translate UniProt ID 'A0A6H2TXZ6' to CARD obtaining all possible translations
if(interactive()){
    getUniProt2CARD('A0A6H2TXZ6', exhaustiveMapping = TRUE, detailedMapping = TRUE)
}
}
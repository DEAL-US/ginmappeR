\name{getCARD2KEGG}
\alias{getCARD2KEGG}
\title{Translates a CARD database ID to KEGG database}
\description{
Retrieves the ID translation of an specific instance of CARD database ID to
the KEGG Genes database.}
\usage{
getCARD2KEGG(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE,
byIdenticalProteins = TRUE, bySimilarGenes = TRUE)
}
\arguments{
  \item{cardId}{valid CARD database ARO ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
  \item{detailedMapping}{logical value (TRUE or FALSE) indicating to return a
  vector of KEGG IDs if FALSE, or a named list if TRUE. This named list details
  the origin of the translated IDs, dividing them in 'DT', if it is a Direct
  Translation, or '1.0', '0.9' or '0.5', if it has been retrieved through any
  of the identical gene clusters. See examples below for further understanding
  }
  \item{byIdenticalProteins}{logical value (TRUE or FALSE) indicating to
  translate also through NCBI Identical Proteins. TRUE value could mean higher
  translation time and FALSE, fewer chances to translate the ID}
  \item{bySimilarGenes}{logical value (TRUE or FALSE) indicating to translate
  also through similar UniProt genes (100\%, 90\% or 50\% of identity). Using
  TRUE value may take a long time, only advised if there is not a direct
  translation possible
  }
}
\value{
Depending on parameter 'detailedMapping': if TRUE it will return a named list
dividing KEGG IDs in 'DT', if it is a Direct Translation, or '1.0', '0.9' or
'0.5', if it has been retrieved through any of the identical gene clusters,
or empty list if no translations were found; if FALSE, returns a vector of
KEGG IDs, empty vector if no translations were found.
}
\examples{
# Translate CARD ID '3000938' to KEGG
if(interactive()){
    getCARD2KEGG('3000938')
}
# Translate CARD ID '3001109' to KEGG
if(interactive()){
    getCARD2KEGG('3001109', detailedMapping = TRUE)}
}

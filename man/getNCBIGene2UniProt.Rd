\name{getNCBIGene2UniProt}
\alias{getNCBIGene2UniProt}
\title{Translates a NCBI Gene database ID to UniProt database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the UniProt database.}
\usage{
getNCBIGene2UniProt(ncbiId, exhaustiveMapping = FALSE,
detailedMapping = FALSE, byIdenticalProteins = TRUE)
}
\arguments{
  \item{ncbiId}{valid NCBI Gene database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
  \item{detailedMapping}{logical value (TRUE or FALSE) indicating whether to
  return a vector of UniProt IDs if FALSE, or a named list if TRUE. This named
  list details the origin of the translated IDs, dividing them in 'DT', if it
  is a Direct Translation, or '1.0' (100\% identity), if it has been retrieved
  through Identical Proteins database. See examples below for further
  understanding}
  \item{byIdenticalProteins}{logical value (TRUE or FALSE) indicating to
  retrieve translations of NCBI Identical Proteins to the given identifier}
}
\value{
Depending on parameter 'detailedMapping': if TRUE it will return a named list
dividing UniProt IDs in 'DT', if it is a Direct Translation, or '1.0', if it
has been retrieved through NCBI Identical Protein database , or empty list if
no translations were found; if FALSE, returns a vector of UniProt IDs, empty
vector if no translations were found.
If an error occurs, returns NULL.
}
\examples{
# Translate NCBI Gene ID '76524190' to UniProt
getNCBIGene2UniProt('76524190')

# Translate NCBI Gene IDs '76524190' and '1272' to UniProt
getNCBIGene2UniProt(c('76524190', '1272'))
}

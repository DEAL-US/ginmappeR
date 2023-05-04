\name{getNCBINucleotide2UniProt}
\alias{getNCBINucleotide2UniProt}
\title{Translates a NCBI Nucleotide database ID to UniProt database}
\description{
Retrieves the ID translation of an specific instance of NCBI Nucleotide database ID to
the UniProt database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBINucleotide2UniProt(ncbiId, exhaustiveMapping = FALSE, byIdenticalProteins = TRUE)
}
\arguments{
  \item{ncbiId}{valid NCBI Nucleotide database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve all possible translations
  or only the first and fastest one}
  \item{byIdenticalProteins}{logical value (TRUE or FALSE) indicating to retrieve translations of
  NCBI Identical Proteins to the given identifier.}
}
\value{
A vector of UniProt IDs, empty if no translations were found.
}
\examples{
# Translate NCBI Nucleotide ID 'AY536519' to UniProt
if(interactive()){
    getNCBINucleotide2UniProt('AY536519')
}
}

\name{getNCBINucleotide2CARD}
\alias{getNCBINucleotide2CARD}
\title{Translates a NCBI Nucleotide database ID to CARD database}
\description{
Retrieves the ID translation of an specific instance of NCBI Nucleotide database
ID to the CARD database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBINucleotide2CARD(ncbiId, exhaustiveMapping = FALSE)
}
\arguments{
  \item{ncbiId}{valid NCBI Nucleotide database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
}
\value{
A vector of CARD IDs, empty if no translations were found.
}
\examples{
# Translate NCBI Nucleotide ID 'AY034138' to CARD
if(interactive()){
    getNCBINucleotide2CARD('AY034138')
}
# Translate NCBI Nucleotide ID 'AJ310778.1' to CARD
if(interactive()){
    getNCBINucleotide2CARD('AJ310778.1')
    }
}

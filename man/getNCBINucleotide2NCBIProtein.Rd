\name{getNCBINucleotide2NCBIProtein}
\alias{getNCBINucleotide2NCBIProtein}
\title{Translates a NCBI Nucleotide database ID to NCBI Protein database}
\description{
Retrieves the ID translation of an specific instance of NCBI Nucleotide database
ID to the NCBI Protein database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBINucleotide2NCBIProtein(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Nucleotide database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
}
\value{
A list of NCBI Protein IDs, empty if no translations were found.
If an error occurs, returns NULL.
}
\examples{
# Translate NCBI Nucleotide ID 'AY536519' to NCBI Protein
getNCBINucleotide2NCBIProtein('AY536519')

# Translate NCBI Nucleotide IDs 'AY536519' and 'X54723.1' to NCBI Protein
getNCBINucleotide2NCBIProtein(c('AY536519', 'X54723.1'))
}

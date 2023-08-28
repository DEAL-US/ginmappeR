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
  \item{id}{valid NCBI Nucleotide database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
}
\value{
A list of NCBI Protein IDs, empty if no translations were found.
}
\examples{
# Translate NCBI Nucleotide ID 'AY536519' to NCBI Protein
if(interactive()){
    getNCBINucleotide2NCBIProtein('AY536519')
}
}

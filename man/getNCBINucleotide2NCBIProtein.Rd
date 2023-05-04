\name{getNCBINucleotide2NCBIProtein}
\alias{getNCBINucleotide2NCBIProtein}
\title{Translates an NCBI Nucleotide database ID to NCBI Protein database}
\description{
Retrieves the ID translation of an specific instance of NCBI Nucleotide database ID to
the NCBI Protein database. It returns a list of the translations or throws an
exception if no translation is found.}
\usage{
getNCBINucleotide2NCBIProtein(id)
}
\arguments{
  \item{id}{valid NCBI Nucleotide database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve all possible translations
  or only the first and fastest one}
}
\value{
A list of translated NCBI Protein IDs.
}
\examples{
# Translate NCBI Nucleotide ID 'AY536519' to NCBI Protein
if(interactive()){
    getNCBINucleotide2NCBIProtein('AY536519')
}
}

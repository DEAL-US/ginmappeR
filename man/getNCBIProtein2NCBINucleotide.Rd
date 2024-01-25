\name{getNCBIProtein2NCBINucleotide}
\alias{getNCBIProtein2NCBINucleotide}
\title{Translates a NCBI Protein database ID to NCBI Nucleotide database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database
ID to the NCBI Nucleotide database. It returns a vector of the translations or
an empty vector if no translation is found. Note that, due to the nature of
NCBI Nucleotide database, a complete genome is often retrieved.}
\usage{
getNCBIProtein2NCBINucleotide(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Protein database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
}
\value{
A list of NCBI Nucleotide IDs, empty if no translations were found.
If an error occurs, returns NULL.
}
\examples{
# Translate NCBI Protein ID 'AAS48620' to NCBI Nucleotide
getNCBIProtein2NCBINucleotide('AAS48620')

# Translate NCBI Protein IDs 'AAS48620' and 'CAA38525' to NCBI Nucleotide
getNCBIProtein2NCBINucleotide(c('AAS48620', 'CAA38525'))
}

\name{getNCBIProtein2NCBIGene}
\alias{getNCBIProtein2NCBIGene}
\title{Translates a NCBI Protein database ID to NCBI Gene database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database
ID to the NCBI Gene database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBIProtein2NCBIGene(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Protein database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
}
\value{
A list of NCBI Gene IDs, empty if no translations were found.
}
\examples{
# Translate NCBI Protein ID 'CAA79696' to NCBI Gene
getNCBIProtein2NCBIGene('CAA79696')
}

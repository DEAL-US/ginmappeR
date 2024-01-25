\name{getNCBIGene2NCBIProtein}
\alias{getNCBIGene2NCBIProtein}
\title{Translates a NCBI Gene database ID to NCBI Protein database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the NCBI Protein database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBIGene2NCBIProtein(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Gene database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
}
\value{
A list of NCBI Protein IDs, empty if no translations were found.
If an error occurs, returns NULL.
}
\examples{
# Translate NCBI Gene ID '76524190' to NCBI Protein
getNCBIGene2NCBIProtein('76524190')

# Translate NCBI Gene IDs '76524190' and '1272' to NCBI Protein
getNCBIGene2NCBIProtein(c('76524190', '1272'))
}

\name{getNCBIGene2NCBIProtein}
\alias{getNCBIGene2NCBIProtein}
\title{Translates an NCBI Gene database ID to NCBI Protein database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the NCBI Protein database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBIGene2NCBIProtein(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Gene database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve all possible translations
  or only the first and fastest one}
}
\value{
A list of NCBI Protein IDs, empty if no translations were found.
}
\examples{
# Translate NCBI Gene ID '76524190' to NCBI Protein
if(interactive()){
    getNCBIGene2NCBIProtein('76524190')
}
}

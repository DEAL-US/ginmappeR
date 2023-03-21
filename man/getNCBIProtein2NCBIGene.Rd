\name{getNCBIProtein2NCBIGene}
\alias{getNCBIProtein2NCBIGene}
\title{Translates an NCBI Protein database ID to NCBI Gene database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database ID to
the NCBI Gene database. It returns a list of the translations or throws an
exception if no translation is found.}
\usage{
getNCBIProtein2NCBIGene(id)
}
\arguments{
  \item{id}{valid NCBI Protein database ID
  }
}
\value{
A list of translated NCBI Gene IDs.
}
\examples{
# Translate NCBI Protein ID 'CAA79696' to NCBI Gene
if(interactive()){
    getNCBIProtein2NCBIGene('CAA79696')
}
}

\name{getNCBIGene2NCBIProtein}
\alias{getNCBIGene2NCBIProtein}
\title{Translates an NCBI Gene database ID to NCBI Protein database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the NCBI Protein database. It returns a list of the translations or throws an
exception if no translation is found.}
\usage{
getNCBIGene2NCBIProtein(id)
}
\arguments{
  \item{id}{valid NCBI Gene database ID
  }
}
\value{
A list of translated NCBI Protein IDs.
}
\examples{
# Translate NCBI Gene ID '76524190' to NCBI Protein
if(interactive()){
    getNCBIGene2NCBIProtein('76524190')
}
}

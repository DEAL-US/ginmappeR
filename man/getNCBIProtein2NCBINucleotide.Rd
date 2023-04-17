\name{getNCBIProtein2NCBINucleotide}
\alias{getNCBIProtein2NCBINucleotide}
\title{Translates an NCBI Protein database ID to NCBI Nucleotide database}
\description{
Retrieves the ID translation of an specific instance of NCBI Protein database ID to
the NCBI Nucleotide database. It returns a list of the translations or throws an
exception if no translation is found.}
\usage{
getNCBIProtein2NCBINucleotide(id)
}
\arguments{
  \item{id}{valid NCBI Protein database ID
  }
}
\value{
A list of translated NCBI Nucleotide IDs.
}
\examples{
# Translate NCBI Protein ID 'CAA79696' to NCBI Nucleotide
if(interactive()){
    getNCBIProtein2NCBINucleotide('CAA79696')
}
}
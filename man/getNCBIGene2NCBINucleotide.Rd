\name{getNCBIGene2NCBINucleotide}
\alias{getNCBIGene2NCBINucleotide}
\title{Translates an NCBI Gene database ID to NCBI Nucleotide database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the NCBI Nucleotide database. It returns a list of the translations or throws an
exception if no translation is found.}
\usage{
getNCBIGene2NCBINucleotide(id)
}
\arguments{
  \item{id}{valid NCBI Gene database ID
  }
}
\value{
A list of translated NCBI Nucleotide IDs.
}
\examples{
# Translate NCBI Gene ID '76524190' to NCBI Nucleotide
if(interactive()){
    getNCBIGene2NCBINucleotide('76524190')
}
}

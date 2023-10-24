\name{getNCBIGene2NCBINucleotide}
\alias{getNCBIGene2NCBINucleotide}
\title{Translates a NCBI Gene database ID to NCBI Nucleotide database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the NCBI Nucleotide database. It returns a vector of the translations or an
empty vector if no translation is found. Note that, due to the nature of
NCBI Nucleotide database, a complete genome is often retrieved.}
\usage{
getNCBIGene2NCBINucleotide(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Gene database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
}
\value{
A vector of NCBI Nucleotide IDs, empty if no translations were found.
}
\examples{
# Translate NCBI Gene ID '76524190' to NCBI Nucleotide
getNCBIGene2NCBINucleotide('76524190')
}

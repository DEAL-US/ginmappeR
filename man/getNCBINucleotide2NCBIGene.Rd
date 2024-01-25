\name{getNCBINucleotide2NCBIGene}
\alias{getNCBINucleotide2NCBIGene}
\title{Translates a NCBI Nucleotide database ID to NCBI Gene database}
\description{
Retrieves the ID translation of an specific instance of NCBI Nucleotide database
ID to the NCBI Gene database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBINucleotide2NCBIGene(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Nucleotide database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
}
\value{
A list of NCBI Gene IDs, empty if no translations were found.
If an error occurs, returns NULL.
}
\examples{
# Translate NCBI Nucleotide ID 'Z21488' to NCBI Gene
getNCBINucleotide2NCBIGene('Z21488')
}

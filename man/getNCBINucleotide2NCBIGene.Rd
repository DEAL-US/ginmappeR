\name{getNCBINucleotide2NCBIGene}
\alias{getNCBINucleotide2NCBIGene}
\title{Translates an NCBI Nucleotide database ID to NCBI Gene database}
\description{
Retrieves the ID translation of an specific instance of NCBI Nucleotide database ID to
the NCBI Gene database. It returns a list of the translations or throws an
exception if no translation is found.}
\usage{
getNCBINucleotide2NCBIGene(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Nucleotide database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve all possible translations
  or only the first and fastest one}
}
\value{
A list of translated NCBI Gene IDs.
}
\examples{
# Translate NCBI Nucleotide ID 'Z21488' to NCBI Gene
if(interactive()){
    getNCBINucleotide2NCBIGene('Z21488')
}
}

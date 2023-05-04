\name{getCARD2NCBINucleotide}
\alias{getCARD2NCBINucleotide}
\title{Translates a CARD database ID to NCBI Nucleotide database}
\description{
Retrieves the ID translation of an specific instance of CARD database ID to
the NCBI Nucleotide database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getCARD2NCBINucleotide(cardId)
}
\arguments{
  \item{cardId}{valid CARD database ARO ID}
}
\value{
A vector of NCBI Nucleotide IDs, empty if no translations were found.
}
\examples{
# Translate CARD ID '3002535' to NCBI Nucleotide
if(interactive()){
    getCARD2NCBINucleotide('3002535')
}
# Translate CARD ID 'ARO:3002535' to NCBI Nucleotide
if(interactive()){
    getCARD2NCBINucleotide('ARO:3002535')
    }
}

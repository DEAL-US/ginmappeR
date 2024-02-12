\name{getCARD2NCBINucleotide}
\alias{getCARD2NCBINucleotide}
\title{Translates a CARD database ID to NCBI Nucleotide database}
\description{
Retrieves the ID translation of an specific instance of CARD database ID to
the NCBI Nucleotide database. Note that, due to the nature of
NCBI Nucleotide database, a complete genome is often retrieved.}
\usage{
getCARD2NCBINucleotide(cardId)
}
\arguments{
  \item{cardId}{valid CARD database ARO ID or a vector of them}
}
\value{
A vector of NCBI Nucleotide IDs.
}
\examples{
# Translate CARD ID '3002535' to NCBI Nucleotide
getCARD2NCBINucleotide('3002535')

# Translate CARD ID 'ARO:3002535' to NCBI Nucleotide
getCARD2NCBINucleotide('ARO:3002535')

# Translate CARD IDs '3002535' and '3003988' to NCBI Nucleotide
getCARD2NCBINucleotide(c('3002535', '3003988'))
}

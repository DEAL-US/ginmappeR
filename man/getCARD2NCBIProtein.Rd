\name{getCARD2NCBIProtein}
\alias{getCARD2NCBIProtein}
\title{Translates a CARD database ID to NCBI Protein database}
\description{
Retrieves the ID translation of an specific instance of CARD database ID to
the NCBI Protein database.}
\usage{
getCARD2NCBIProtein(cardId)
}
\arguments{
  \item{cardId}{valid CARD database ARO ID or a vector of them}
}
\value{
A character vector of the same length as \code{cardId}.
}
\examples{
# Translate CARD ID '3002535' to NCBI Protein
getCARD2NCBIProtein('3002535')

# Translate CARD ID 'ARO:3002535' to NCBI Protein
getCARD2NCBIProtein('ARO:3002535')

# Translate CARD IDs '3002535' and '3003988' to NCBI Protein
getCARD2NCBIProtein(c('3002535', '3003988'))
}

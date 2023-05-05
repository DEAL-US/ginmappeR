\name{getCARD2NCBIProtein}
\alias{getCARD2NCBIProtein}
\title{Translates a CARD database ID to NCBI Protein database}
\description{
Retrieves the ID translation of an specific instance of CARD database ID to
the NCBI Protein database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getCARD2NCBIProtein(cardId)
}
\arguments{
  \item{cardId}{valid CARD database ARO ID}
}
\value{
A vector of NCBI Protein IDs, empty if no translations were found.
}
\examples{
# Translate CARD ID '3002535' to NCBI Protein
if(interactive()){
    getCARD2NCBIProtein('3002535')
}
# Translate CARD ID 'ARO:3002535' to NCBI Protein
if(interactive()){
    getCARD2NCBIProtein('ARO:3002535')
    }
}
\name{getCARD2NCBIGene}
\alias{getCARD2NCBIGene}
\title{Translates a CARD database ID to NCBI Gene database}
\description{
Retrieves the ID translation of an specific instance of CARD database ID to
the NCBI Gene database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getCARD2NCBIGene(cardId, exhaustiveMapping = FALSE)
}
\arguments{
  \item{cardId}{valid CARD database ARO ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve all possible translations
  or only the first and fastest one}
}
\value{
A vector of NCBI Gene IDs, empty if no translations were found.
}
\examples{
# Translate CARD ID '3002524' to NCBI Gene
if(interactive()){
    getCARD2NCBIGene('3002524')
}
# Translate CARD ID 'ARO:3002524' to NCBI Gene
if(interactive()){
    getCARD2NCBIGene('ARO:3002524', TRUE)
    }
}

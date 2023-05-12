\name{getCARD2UniProt}
\alias{getCARD2UniProt}
\title{Translates a CARD database ID to UniProt database}
\description{
Retrieves the ID translation of an specific instance of CARD database ID to
the UniProt database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getCARD2UniProt(cardId, exhaustiveMapping = FALSE)
}
\arguments{
  \item{cardId}{valid CARD database ARO ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
}
\value{
A vector of UniProt IDs, empty if no translations were found.
}
\examples{
# Translate CARD ID '3002867' to UniProt
if(interactive()){
    getCARD2UniProt('3002867')
}
}

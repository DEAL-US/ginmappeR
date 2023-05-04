\name{getNCBIGene2CARD}
\alias{getNCBIGene2CARD}
\title{Translates a NCBI Gene database ID to CARD database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the CARD database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getNCBIGene2CARD(ncbiId)
}
\arguments{
  \item{ncbiId}{valid NCBI Gene database ID}
}
\value{
A vector of CARD IDs, empty if no translations were found.
}
\examples{
# Translate NCBI Gene ID '3510143' to CARD
if(interactive()){
    getNCBIGene2CARD('3510143')
}
}

\name{getNCBIGene2NCBIProtein}
\alias{getNCBIGene2NCBIProtein}
\title{Translates a NCBI Gene database ID to NCBI Protein database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the NCBI Protein database.}
\usage{
getNCBIGene2NCBIProtein(id, exhaustiveMapping=FALSE)
}
\arguments{
  \item{id}{valid NCBI Gene database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (\code{TRUE} or \code{FALSE})
  indicating to retrieve all possible translations or only the first and
  fastest one.}
}
\value{
A character vector of the same length as \code{id}. If parameter
\code{exhaustiveMapping} is \code{TRUE}, it will return a list instead of a
character vector to avoid mixing the identifiers.
}
\examples{
\donttest{
# Translate NCBI Gene ID '76524190' to NCBI Protein
getNCBIGene2NCBIProtein('76524190')

# Translate NCBI Gene IDs '76524190' and '1272' to NCBI Protein
getNCBIGene2NCBIProtein(c('76524190', '1272'))

# Get all possible translations of IDs '76524190', '1272'
getNCBIGene2NCBIProtein(c('76524190', '1272'), exhaustiveMapping = TRUE)
}
}

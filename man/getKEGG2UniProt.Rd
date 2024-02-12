\name{getKEGG2UniProt}
\alias{getKEGG2UniProt}
\title{Translates a KEGG Genes database ID to UniProt database}
\description{
Retrieves the ID translation of an specific instance of KEGG Genes database ID
to the UniProt database.}
\usage{
getKEGG2UniProt(keggId, exhaustiveMapping = FALSE)
}
\arguments{
  \item{keggId}{valid KEGG Genes database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (\code{TRUE} or \code{FALSE})
  indicating to retrieve all possible translations or only the first and
  fastest one.}
}
\value{
A character vector of the same length as \code{keggId}. If parameter
\code{exhaustiveMapping} is \code{TRUE}, it will return a list instead of a
character vector to avoid mixing the identifiers.
}
\examples{
\donttest{
# Translate KEGG Genes ID 'llo:LLO_2673' to UniProt
getKEGG2UniProt('llo:LLO_2673')

# Translate KEGG IDs 'llo:LLO_2673' and 'abh:M3Q_2831' to NCBI UniProt
getKEGG2UniProt(c('llo:LLO_2673', 'abh:M3Q_2831'))

# Get all possible translations of IDs 'llo:LLO_2673', 'abh:M3Q_2831'
getKEGG2UniProt(c('llo:LLO_2673', 'abh:M3Q_2831'), exhaustiveMapping = TRUE)
}
}

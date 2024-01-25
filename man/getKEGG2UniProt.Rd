\name{getKEGG2UniProt}
\alias{getKEGG2UniProt}
\title{Translates a KEGG Genes database ID to UniProt database}
\description{
Retrieves the ID translation of an specific instance of KEGG Genes database ID
to the UniProt database. It returns a vector of the translations or an
empty vector if no translation is found.}
\usage{
getKEGG2UniProt(keggId)
}
\arguments{
  \item{keggId}{valid KEGG Genes database ID or a vector of them}
}
\value{
A vector of UniProt IDs, empty if no translations were found.
If an error occurs, returns NULL.
}
\examples{
# Translate KEGG Genes ID 'llo:LLO_2673' to UniProt
getKEGG2UniProt('llo:LLO_2673')

# Translate KEGG IDs 'llo:LLO_2673' and 'abh:M3Q_2831' to NCBI UniProt
getKEGG2UniProt(c('llo:LLO_2673', 'abh:M3Q_2831'))
}

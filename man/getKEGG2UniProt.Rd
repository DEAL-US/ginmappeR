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
  \item{keggId}{valid KEGG Genes database ID}
}
\value{
A vector of UniProt IDs, empty if no translations were found.
}
\examples{
# Translate KEGG Genes ID 'lacy:A4V08_06425' to UniProt
if(interactive()){
    getKEGG2UniProt('lacy:A4V08_06425')
}
# Translate KEGG Genes ID 'aag:5579347' to UniProt
if(interactive()){
    getKEGG2UniProt('aag:5579347')
}
}

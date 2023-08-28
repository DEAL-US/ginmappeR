\name{getKEGG2NCBIProtein}
\alias{getKEGG2NCBIProtein}
\title{Translates a KEGG database ID to NCBI Protein database}
\description{
Retrieves the ID translation of an specific instance of KEGG Genes database ID
to the NCBI Protein database.}
\usage{
getKEGG2NCBIProtein(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE,
bySimilarGenes = TRUE)
}
\arguments{
  \item{keggId}{valid KEGG Genes database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
  \item{detailedMapping}{logical value (TRUE or FALSE) indicating to return a
  vector of NCBI Protein IDs if FALSE, or a named list if TRUE. This named list
  details the origin of the translated IDs, dividing them in 'DT', if it is a
  Direct Translation, or '1.0', '0.9' or '0.5', if it has been retrieved through
  any of the identical gene clusters.
  }
  \item{bySimilarGenes}{logical value (TRUE or FALSE) indicating to translate
  also through similar UniProt genes (100\%, 90\% or 50\% of identity). Using
  TRUE value may take a long time, only advised if there is not a direct
  translation possible
  }
}
\value{
Depending on parameter 'detailedMapping': if TRUE it will return a named list
dividing NCBI Protein IDs in 'DT', if it is a Direct Translation, or '1.0',
'0.9' or '0.5', if it has been retrieved through any of the identical gene
clusters, or empty list if no translations were found; if FALSE, returns a
vector of NCBI Protein IDs, empty vector if no translations were found.
}
\examples{
# Translate KEGG ID 'llo:LLO_2673' to NCBI Protein
if(interactive()){
getKEGG2NCBIProtein('llo:LLO_2673')
}
}

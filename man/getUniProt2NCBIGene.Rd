\name{getUniProt2NCBIGene}
\alias{getUniProt2NCBIGene}
\title{Translates an UniProt database ID to NCBI Gene database}
\description{
Retrieves the ID translation of an specific instance of UniProt database ID to
the NCBI Gene database.}
\usage{
getUniProt2NCBIGene(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE,
bySimilarGenes = TRUE)
}
\arguments{
  \item{upId}{valid UniProt database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
  \item{detailedMapping}{logical value (TRUE or FALSE) indicating to return a
  vector of NCBI Gene IDs if FALSE, or a named list if TRUE. This named list
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
dividing NCBI Gene IDs in 'DT', if it is a Direct Translation, or '1.0', '0.9'
or '0.5', if it has been retrieved through any of the identical gene clusters,
or empty list if no translations were found; if FALSE, returns a vector of NCBI
Gene IDs, empty vector if no translations were found.
}
\examples{
# Translate UniProt ID 'A0SNL9' to NCBI Gene
if(interactive()){
getUniProt2NCBIGene('A0SNL9')
}
}

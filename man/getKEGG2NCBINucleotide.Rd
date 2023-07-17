\name{getKEGG2NCBINucleotide}
\alias{getKEGG2NCBINucleotide}
\title{Translates a KEGG database ID to NCBI Nucleotide database}
\description{
Retrieves the ID translation of an specific instance of KEGG Genes database ID
to the NCBI Nucleotide database. Note that, due to the nature of
NCBI Nucleotide database, a complete genome is often retrieved.}
\usage{
getKEGG2NCBINucleotide(keggId, exhaustiveMapping = FALSE,
detailedMapping = FALSE, bySimilarGenes = TRUE)
}
\arguments{
  \item{keggId}{valid KEGG Genes database ID}
  \item{exhaustiveMapping}{logical value (TRUE or FALSE) indicating to retrieve
  all possible translations or only the first and fastest one}
  \item{detailedMapping}{logical value (TRUE or FALSE) indicating to return a
  vector of NCBI Nucleotide IDs if FALSE, or a named list if TRUE. This named
  list details the origin of the translated IDs, dividing them in 'DT', if it
  is a Direct Translation, or '1.0', '0.9' or '0.5', if it has been retrieved
  through any of the identical gene clusters.
  }
  \item{bySimilarGenes}{logical value (TRUE or FALSE) indicating to translate
  also through similar UniProt genes (100\%, 90\% or 50\% of identity). Using
  TRUE value may take a long time, only advised if there is not a direct
  translation possible
  }
}
\value{
Depending on parameter 'detailedMapping': if TRUE it will return a named list
dividing NCBI Nucleotide IDs in 'DT', if it is a Direct Translation, or '1.0',
'0.9' or '0.5', if it has been retrieved through any of the identical gene
clusters, or empty list if no translations were found; if FALSE, returns a
vector of NCBI Nucleotide IDs, empty vector if no translations were found.
}
\examples{
# Translate KEGG ID 'llo:LLO_2673' to NCBI Nucleotide
if(interactive()){
getKEGG2NCBINucleotide('llo:LLO_2673')
}
# Translate KEGG ID 'llo:LLO_2673' to NCBI Nucleotide obtaining all possible
# translations
if(interactive()){
getKEGG2NCBINucleotide('llo:LLO_2673', exhaustiveMapping = TRUE,
detailedMapping = TRUE)
}
}

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
  \item{keggId}{valid KEGG Genes database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (\code{TRUE} or \code{FALSE})
  indicating to retrieve all possible translations or only the first and
  fastest one.}
  \item{detailedMapping}{logical value (\code{TRUE} or \code{FALSE}) indicating
  to return a vector of NCBI Nucleotide IDs if \code{FALSE}, or a named vector
  if \code{TRUE}. This named vector details the origin of the translated IDs,
  dividing them in \code{DT}, if it is a Direct Translation, or \code{1.0},
  \code{0.9} or \code{0.5}, if it has been retrieved through any of the
  identical gene clusters.
  }
  \item{bySimilarGenes}{logical value (\code{TRUE} or \code{FALSE}) indicating
  to translate also through similar UniProt genes (100\%, 90\% or 50\% of
  identity). Using \code{TRUE} value may take a long time, only advised if
  there is not a direct translation possible.
  }
}
\value{
A character vector of the same length as \code{keggId}. If parameter
\code{exhaustiveMapping} is \code{TRUE}, it will return a list instead of a
character vector to avoid mixing the identifiers.
}
\examples{
\donttest{
# Translate KEGG ID 'llo:LLO_2673' to NCBI Nucleotide
getKEGG2NCBINucleotide('llo:LLO_2673')

# Translate KEGG IDs 'llo:LLO_2673' and 'abh:M3Q_2831' to NCBI Nucleotide
getKEGG2NCBINucleotide(c('llo:LLO_2673', 'abh:M3Q_2831'))

# Get all possible translations of IDs 'aag:5579347', 'abh:M3Q_2831' with
# detailed origin
getKEGG2NCBINucleotide(c('aag:5579347', 'abh:M3Q_2831'),
exhaustiveMapping = TRUE, detailedMapping = TRUE)
}
}

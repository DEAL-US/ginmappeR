\name{getKEGG2NCBIGene}
\alias{getKEGG2NCBIGene}
\title{Translates a KEGG database ID to NCBI Gene database}
\description{
Retrieves the ID translation of an specific instance of KEGG Genes database ID
to the NCBI Gene database.}
\usage{
getKEGG2NCBIGene(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE,
bySimilarGenes = TRUE)
}
\arguments{
  \item{keggId}{valid KEGG Genes database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (\code{TRUE} or \code{FALSE})
  indicating to retrieve all possible translations or only the first and
  fastest one.}
  \item{detailedMapping}{logical value (\code{TRUE} or \code{FALSE}) indicating
  to return a vector of NCBI Gene IDs if \code{FALSE}, or a named vector if
  \code{TRUE}. This named vector details the origin of the translated IDs,
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
# Translate KEGG ID 'aag:5579347' to NCBI Gene
getKEGG2NCBIGene('aag:5579347', detailedMapping = TRUE)
\donttest{
# Translate KEGG IDs 'aag:5579347' and 'abh:M3Q_2831' to NCBI Gene
getKEGG2NCBIGene(c('aag:5579347', 'abh:M3Q_2831'))

# Get all possible translations of IDs 'aag:5579347', 'abh:M3Q_2831' with
# detailed origin
getKEGG2NCBIGene(c('aag:5579347', 'abh:M3Q_2831'), exhaustiveMapping = TRUE,
detailedMapping = TRUE)
}
}

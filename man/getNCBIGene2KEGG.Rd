\name{getNCBIGene2KEGG}
\alias{getNCBIGene2KEGG}
\title{Translates a NCBI Gene database ID to KEGG database}
\description{
Retrieves the ID translation of an specific instance of NCBI Gene database ID to
the KEGG Genes database.}
\usage{
getNCBIGene2KEGG(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE,
byIdenticalProteins = TRUE, bySimilarGenes = TRUE)
}
\arguments{
  \item{ncbiId}{valid NCBI Gene database ID or a vector of them}
  \item{exhaustiveMapping}{logical value (\code{TRUE} or \code{FALSE})
  indicating to retrieve all possible translations or only the first and
  fastest one.}
  \item{detailedMapping}{logical value (\code{TRUE} or \code{FALSE}) indicating
  to return a vector of KEGG IDs if \code{FALSE}, or a named vector if
  \code{TRUE}. This named vector details the origin of the translated IDs,
  dividing them in \code{DT}, if it is a Direct Translation, or \code{1.0},
  \code{0.9} or \code{0.5}, if it has been retrieved through any of the
  identical gene clusters.
  }
  \item{byIdenticalProteins}{logical value (\code{TRUE} or \code{FALSE})
  indicating to translate also through NCBI Identical Proteins. \code{TRUE}
  value could mean higher translation time and \code{FALSE}, fewer chances to
  translate the ID}
  \item{bySimilarGenes}{logical value (\code{TRUE} or \code{FALSE}) indicating
  to translate also through similar UniProt genes (100\%, 90\% or 50\% of
  identity). Using \code{TRUE} value may take a long time, only advised if
  there is not a direct translation possible.
  }
}
\value{
A character vector of the same length as \code{ncbiId}. If parameter
\code{exhaustiveMapping} is \code{TRUE}, it will return a list instead of a
character vector to avoid mixing the identifiers.
}
\examples{
# Translate NCBI Gene ID '76524190' to KEGG
getNCBIGene2KEGG('AFU38919', bySimilarGenes = FALSE, byIdenticalProteins = FALSE)
\donttest{
# Translate NCBI Gene IDs '76524190' and 'AFU38919' to KEGG
getNCBIGene2KEGG(c('76524190', 'AFU38919'))

# Get all possible translations of IDs '76524190' and 'AFU38919' with
# detailed origin
getNCBIGene2KEGG(c('76524190', 'AFU38919'), detailedMapping = TRUE,
exhaustiveMapping = TRUE)
}
}

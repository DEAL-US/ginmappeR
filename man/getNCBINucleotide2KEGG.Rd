\name{getNCBINucleotide2KEGG}
\alias{getNCBINucleotide2KEGG}
\title{Translates a NCBI Nucleotide database ID to KEGG database}
\description{
Retrieves the ID translation of an specific instance of NCBI Nucleotide database
ID to the KEGG Genes database.}
\usage{
getNCBINucleotide2KEGG(ncbiId, exhaustiveMapping = FALSE,
detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE)
}
\arguments{
  \item{ncbiId}{valid NCBI Nucleotide database ID or a vector of them}
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
\donttest{
# Translate NCBI Nucleotide ID 'AY536519' to KEGG
getNCBINucleotide2KEGG('AY536519')

# Translate NCBI Nucleotide IDs 'AY536519' and 'NC_013861.1' to KEGG
getNCBINucleotide2KEGG(c('AY536519', 'NC_013861.1'))

# Get all possible translations of IDs 'AY536519', 'NC_013861.1'
# with detailed origin
getNCBINucleotide2KEGG(c('AY536519', 'NC_013861.1'), detailedMapping = TRUE,
exhaustiveMapping = TRUE)
}
}

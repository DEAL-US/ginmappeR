\name{getNCBIIdenticalProteins}
\alias{getNCBIIdenticalProteins}
\title{Search NCBI Identical Protein database by a given NCBI ID}
\description{
Retrieves identical NCBI proteins to that of a given ID.}
\usage{
getNCBIIdenticalProteins(ncbiId, format='ids')
}
\arguments{
    \item{ncbiId}{valid NCBI ID or a vector of them}
    \item{format}{return format of the response. Default value \code{ids}
    returns a vector of the identical proteins IDs; \code{dataframe} value
    returns a dataframe with different information about the found identical
    proteins}
}
\value{
Depending on parameter \code{format}, it will return a list of same legth as
\code{ncbiId} with identical NCBI Protein IDs or with dataframes containing
different information about the identical proteins found: \code{Id},
\code{Source}, \code{Nucleotide.Accession}, \code{Start}, \code{Stop},
\code{Strand}, \code{Protein}, \code{Protein.Name}, \code{Organism},
\code{Strain} and \code{Assembly}.
}
\examples{
# Get IDs on identical proteins to the one with ID 'AHA80958'
getNCBIIdenticalProteins('AHA80958', format = 'ids')
\donttest{
# Get IDs on identical proteins to the ones with IDs 'AHA80958' and 'CAA76794'
getNCBIIdenticalProteins(c('AHA80958', 'CAA76794'), format = 'ids')
}
}

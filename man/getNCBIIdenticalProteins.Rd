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
    \item{format}{return format of the response. Default value 'ids' returns
    a vector of the identical proteins IDs; 'dataframe' value returns a
    dataframe with different information about the found identical proteins}

}
\value{
Depending on parameter 'format', it will return a vector of identical NCBI
Protein IDs or a dataframe containing different information about the identical
proteins found: 'Id', 'Source', 'Nucleotide.Accession', 'Start', 'Stop',
'Strand', 'Protein', 'Protein.Name', 'Organism', 'Strain' and 'Assembly'.
}
\examples{
# Get IDs on identical proteins to the one with ID 'AHA80958'
getNCBIIdenticalProteins('AHA80958', format = 'ids')
}

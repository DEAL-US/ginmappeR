\name{getNCBIIdenticalProteins}
\alias{getNCBIIdenticalProteins}
\title{Searches NCBI Identical Protein Groups database by a given NCBI ID}
\description{
Retrieves a list of identical proteins to that of a given ID. It throws an
exception if no identical proteins are found.}
\usage{
getNCBIIdenticalProteins(ncbiId, format='ids')
}
\arguments{
    \item{ncbiId}{valid NCBI ID}
    \item{format}{return format of the response. Default value 'ids' returns
    a list of the identical proteins IDs; 'dataframe' value returns a dataframe
    with different information about the found identical proteins}

}
\value{
Depending on parameter 'format', it will return a list of identical NCBI Protein IDs
or a dataframe containing different information about the identical proteins found:
'Id', 'Source', 'Nucleotide.Accession', 'Start', 'Stop', 'Strand', 'Protein', 'Protein.Name',
'Organism', 'Strain' and 'Assembly'.
}
\examples{
# Get IDs on identical proteins to the one with ID 'AHA80958'
if(interactive()){
    getNCBIIdenticalProteins('AHA80958', format = 'ids')
}
# Get data on identical proteins to the one with ID 'AHA80958'
if(interactive()){
    getNCBIIdenticalProteins('AHA80958', format = 'dataframe')

}
}

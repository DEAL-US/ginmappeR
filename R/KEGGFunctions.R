############################
# KEGG database to UniProt #
############################

getKEGG2UniProt <- function(keggId){
    .checkKEGGIdExists(keggId)

    query <- keggConv('uniprot', keggId)
    if(identical(query, character(0))|identical(query, "")){
        return(character(0))
    }else{
        splittedQuery <- strsplit(query,':')
        return(unname(unlist(splittedQuery)[2*(1:length(query))]))
    }
}

###################################
# KEGG database to NCBI databases #
###################################

# 1st: direct translation through keggConv
# 2st: translate first to UniProt and use
# getUniProt2NCBIxxx methods to get to NCBIxxx
# To consider: must check what NCBI database
# have to translate and return its identifier.

.getKEGG2NCBIDT <- function(keggId){

    query <- keggConv('ncbi-geneid', keggId)
    query <- c(query, keggConv('ncbi-proteinid', keggId))

    if(identical(query, character(0))|identical(query, "")){
        return(character(0))
    }else{
        splittedQuery <- strsplit(query,':')
        return(unname(unlist(splittedQuery)[2*(1:length(query))]))
    }
}

.getKEGG2NCBITUP <- function(keggId, ncbiDB, exhaustiveMapping = FALSE, bySimilarGenes = TRUE){

    result <- list()
    upIDs <- getKEGG2UniProt(keggId)

    if(length(upIDs)>0){

        for(upId in upIDs){
            ncbiIDs <- switch (ncbiDB,
                 'protein' = getUniProt2NCBIProtein(upId, exhaustiveMapping = exhaustiveMapping, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes),
                 'nucleotide' = getUniProt2NCBINucleotide(upId, exhaustiveMapping = exhaustiveMapping, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes),
                 'gene' = getUniProt2NCBIGene(upId, exhaustiveMapping = exhaustiveMapping, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes),
            )

            if(length(ncbiIDs)>0){
                result <- .mergeNamedLists(result, ncbiIDs)
            }
        }
    }
    return(lapply(result, unique))
}



# TODO: Warnear: using bySimilarGenes=TRUE with exhaustiveMapping=TRUE takes a long time
# .getKEGG2NCBI <- function(){
#
# }






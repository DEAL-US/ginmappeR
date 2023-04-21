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

# Direct translation through keggConv
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

# Through UniProt method: translate first to UniProt and use
# getUniProt2NCBIxxx methods to get to NCBIxxx
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

# TODO: Warning to be included in documentation: using bySimilarGenes=TRUE with exhaustiveMapping=TRUE takes a long time
.getKEGG2NCBI <- function(keggId, ncbiDB, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkKEGGIdExists(keggId)
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')

    result <- list()
    translationsDT <- .getKEGG2NCBIDT(keggId)
    translationsDT <- translationsDT[which(sapply(translationsDT, .checkIdInNCBIDataBase, ncbiDB=ncbiDB))]

    # First, try to retrieve a direct translation
    if(!identical(translationsDT, character(0))){
        if(!exhaustiveMapping){
            if(!detailedMapping){
                return(translationsDT[1])
            }else{
                result[['DT']] <- translationsDT[1]
            }
            return(result)
        }else{
            result[['DT']] <- translationsDT
        }
    }

    # Then, try to translate through UniProt database
    result <- .mergeNamedLists(result, .getKEGG2NCBITUP(keggId, ncbiDB, exhaustiveMapping = exhaustiveMapping,
                                    bySimilarGenes = bySimilarGenes))

    if(!detailedMapping){
        result <- unique(unlist(result, recursive = FALSE, use.names = FALSE))
        if(is.null(result)){result <- character(0)}
    }else{
        result <- lapply(result, unique)
    }
    return(result)
}

getKEGG2NCBIProtein <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    return(.getKEGG2NCBI(keggId, 'protein', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                            bySimilarGenes = bySimilarGenes))
}

getKEGG2NCBINucleotide <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    return(.getKEGG2NCBI(keggId, 'nucleotide', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                            bySimilarGenes = bySimilarGenes))
}

getKEGG2NCBIGene <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    return(.getKEGG2NCBI(keggId, 'gene', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                            bySimilarGenes = bySimilarGenes))
}


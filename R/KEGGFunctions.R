############################
# KEGG database to UniProt #
############################

getKEGG2UniProt <- function(keggId){

    tryCatch(
        {
            .checkKEGGIdExists(keggId)

            query <- keggConv('uniprot', keggId)
            if(identical(query, character(0))|identical(query, "")){
                return(character(0))
            }else{
                splittedQuery <- strsplit(query,':')
                return(unname(unlist(splittedQuery)[2*(1:length(query))]))
            }
        },
        error = function(e) {
            warning(conditionMessage(e), "\n", call. = FALSE, noBreaks. = TRUE, immediate. = TRUE)
            return(NULL)
        }
    )
}

###################################
# KEGG database to NCBI databases #
###################################

# Direct translation through keggConv
.getKEGG2NCBIDT <- function(keggId){

    query <- keggConv('ncbi-geneid', keggId)
    query <- c(query, keggConv('ncbi-proteinid', keggId))

    splittedQuery <- strsplit(query,':')
    splittedQuery <- unname(unlist(splittedQuery)[2*(1:length(query))])

    aux <- keggGet(keggId)
    if('DBLINKS' %in% names(aux[[1]])){
        aux <- unlist(lapply(strsplit(aux[[1]]$DBLINKS[grep("^NCBI", aux[[1]]$DBLINKS)], ":"), function(x) trimws(x[2])))
        if(is.null(aux)){aux <- character(0)}
        splittedQuery <- unique(c(splittedQuery, aux))
    }

    if(identical(splittedQuery, character(0))|is.null(splittedQuery)|identical(splittedQuery, "")){
        return(character(0))
    }else{
        return(splittedQuery)
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

.getKEGG2NCBI <- function(keggId, ncbiDB, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkKEGGIdExists(keggId)

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
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    tryCatch(
        {
            return(.getKEGG2NCBI(keggId, 'protein', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                                 bySimilarGenes = bySimilarGenes))
        },
        error = function(e) {
            warning(conditionMessage(e), "\n", call. = FALSE, noBreaks. = TRUE, immediate. = TRUE)
            return(NULL)
        }
    )
}

getKEGG2NCBINucleotide <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    tryCatch(
        {
            return(.getKEGG2NCBI(keggId, 'nucleotide', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                                 bySimilarGenes = bySimilarGenes))
        },
        error = function(e) {
            warning(conditionMessage(e), "\n", call. = FALSE, noBreaks. = TRUE, immediate. = TRUE)
            return(NULL)
        }
    )
}

getKEGG2NCBIGene <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    tryCatch(
        {
            return(.getKEGG2NCBI(keggId, 'gene', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                                 bySimilarGenes = bySimilarGenes))
        },
        error = function(e) {
            warning(conditionMessage(e), "\n", call. = FALSE, noBreaks. = TRUE, immediate. = TRUE)
            return(NULL)
        }
    )
}

#########################
# KEGG database to CARD #
#########################

# Auxiliar function to retrieve intermediate NCBI ids and check their translations
# to CARD on the go
.getKEGG2CARDexhaustiveMappingAux <- function(keggId, ncbiDB, bySimilarGenes){
    aroIndex <- read.csv(paste(options("cardPath")[[1]],'/card-data/aro_index.tsv',sep=''), sep='\t')

    # This function translates a bunch of NCBI ids to CARD
    .auxNCBI2CARD <- function(translations, ncbiDB){
        if(identical(ncbiDB, 'protein')){
            translations <- lapply(translations, FUN = function(x) unlist(sapply(x, USE.NAMES = FALSE, FUN = function(auxId){
                auxId <- strsplit(auxId,'.', fixed = TRUE)[[1]][[1]]
                return(aroIndex[gsub("\\..*", "", aroIndex$Protein.Accession) == auxId, "ARO.Accession"])
            })))
        }else{
            translations <- lapply(translations, FUN = function(x) unlist(sapply(x, USE.NAMES = FALSE, FUN = function(auxId){
                auxId <- strsplit(auxId,'.', fixed = TRUE)[[1]][[1]]
                return(aroIndex[gsub("\\..*", "", aroIndex$DNA.Accession) == auxId, "ARO.Accession"])
            })))
        }
        return(translations)
    }

    # Direct translation
    translationsDT <- .getKEGG2NCBIDT(keggId)
    if(!identical(translationsDT, character(0))){
        translationsDT <- lapply(.auxNCBI2CARD(list('DT' = translationsDT),ncbiDB), unique)
        translationsDT <- translationsDT[lengths(translationsDT) > 0]
        if(length(translationsDT[['DT']]) > 0){
            return(translationsDT)
        }
    }

    # Through UniProt
    translationsTUP <- .getKEGG2NCBITUP(keggId, ncbiDB, exhaustiveMapping = TRUE, bySimilarGenes = bySimilarGenes)
    if(length(translationsTUP)>0){
        translationsTUP <- lapply(.auxNCBI2CARD(translationsTUP, ncbiDB), unique)
        translationsTUP <- translationsTUP[lengths(translationsTUP) > 0]
        if(length(translationsTUP)>0){
            return(translationsTUP)
        }else{
            return(list())
        }
    }
}


# Added bySimilarGenes param to deactivate it in case it takes long time to finish, or activate it if no
# translation is found
getKEGG2CARD <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkIfCARDIsDownloaded()
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')

    tryCatch(
        {
            .checkKEGGIdExists(keggId)

            aroIndex <- read.csv(paste(options("cardPath")[[1]],'/card-data/aro_index.tsv',sep=''), sep='\t')
            result <- list()

            # Handle not exhaustiveMapping case, trying to get the fastest result present in CARD
            if(!exhaustiveMapping){
                proteinId <-.getKEGG2CARDexhaustiveMappingAux(keggId, 'protein', bySimilarGenes)
                if(length(proteinId)>0){
                    proteinId <- Filter(Negate(is.null), proteinId[c("DT", "1.0", "0.9", "0.5")])
                    result[names(proteinId)[1]] <- proteinId[[1]][[1]]
                }else{
                    nucleotideId <-.getKEGG2CARDexhaustiveMappingAux(keggId, 'nucleotide', bySimilarGenes)
                    if(length(nucleotideId)>0){
                        nucleotideId <- Filter(Negate(is.null), nucleotideId[c("DT", "1.0", "0.9", "0.5")])
                        result[names(nucleotideId)[1]] <- nucleotideId[[1]][[1]]
                    }
                }
            }else{
                # Handle exhaustiveMapping case, trying to retrieve all possible cases
                proteinId <- getKEGG2NCBIProtein(keggId, exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes)
                if(length(proteinId)>0){
                    result <- lapply(proteinId, FUN = function(x) unlist(sapply(x, USE.NAMES = FALSE, FUN = function(auxId){
                        auxId <- strsplit(auxId,'.', fixed = TRUE)[[1]][[1]]
                        return(aroIndex[gsub("\\..*", "", aroIndex$Protein.Accession) == auxId, "ARO.Accession"])
                    })))
                    result <- result[lengths(result) > 0]
                }

                nucleotideId <- getKEGG2NCBINucleotide(keggId, exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes)
                if(length(nucleotideId)>0){
                    nucleotideId <- lapply(nucleotideId, FUN = function(x) unlist(sapply(x, USE.NAMES = FALSE, FUN = function(auxId){
                        auxId <- strsplit(auxId,'.', fixed = TRUE)[[1]][[1]]
                        return(aroIndex[gsub("\\..*", "", aroIndex$DNA.Accession) == auxId, "ARO.Accession"])
                    })))
                }
                nucleotideId <- nucleotideId[lengths(nucleotideId) > 0]
                result <- .mergeNamedLists(result, nucleotideId)
            }

            # Finally, parse the output according to detailedMapping parameter
            if(!detailedMapping){
                result <- unique(unlist(result, recursive = FALSE, use.names = FALSE))
                if(is.null(result)){result <- character(0)}
            }else{
                result <- lapply(result, unique)
                result <- result[lengths(result) > 0]
                if(length(result)==0){result <- list()}
            }
            return(result)
        },
        error = function(e) {
            warning(conditionMessage(e), "\n", call. = FALSE, noBreaks. = TRUE, immediate. = TRUE)
            return(NULL)
        }
    )
}


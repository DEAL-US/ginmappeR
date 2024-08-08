############################
# KEGG database to UniProt #
############################

.getKEGG2UniProt <- function(keggId, exhaustiveMapping = FALSE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .updateProgressBar("Connecting to KEGG web services")
    .updateProgressBar(paste("Connecting to KEGG API and translating id ", keggId," to Uniprot", sep=""))

    tryCatch(
        {
            .checkKEGGIdExists(keggId)

            query <- keggConv('uniprot', keggId)
            if(identical(query, character(0)) || identical(query, "")){
                return(character(0))
            }else{
                splittedQuery <- strsplit(query,':')
                splittedQuery <- unname(unlist(splittedQuery)[2*(1:length(query))])
                if(!exhaustiveMapping){
                    splittedQuery <- splittedQuery[[1]]
                }
                return(splittedQuery)
            }
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
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

    if(identical(splittedQuery, character(0)) || is.null(splittedQuery) || identical(splittedQuery, "")){
        return(character(0))
    }else{
        return(splittedQuery)
    }
}


# Through UniProt method: translate first to UniProt and use
# getUniProt2NCBIxxx methods to get to NCBIxxx
.getKEGG2NCBITUP <- function(keggId, ncbiDB, exhaustiveMapping = FALSE, bySimilarGenes = TRUE){

    result <- list()
    upIDs <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getKEGG2UniProt(keggId)))

    if(length(upIDs)>0){

        for(upId in upIDs){
            ncbiIDs <- switch (ncbiDB,
                 'protein' = .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getUniProt2NCBIProtein(upId, exhaustiveMapping = exhaustiveMapping, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes))),
                 'nucleotide' = .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getUniProt2NCBINucleotide(upId, exhaustiveMapping = exhaustiveMapping, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes))),
                 'gene' = .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getUniProt2NCBIGene(upId, exhaustiveMapping = exhaustiveMapping, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes)))
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

    .updateProgressBar(paste("Trying to get a direct translation for KEGG id ", keggId," to NCBI", sep=""))

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

    .updateProgressBar(paste("Translating KEGG id ", keggId," to NCBI through similar genes strategy", sep=""))

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

.getKEGG2NCBIProtein <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    tryCatch(
        {
            return(.getKEGG2NCBI(keggId, 'protein', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                                 bySimilarGenes = bySimilarGenes))
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}

.getKEGG2NCBINucleotide <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    tryCatch(
        {
            return(.getKEGG2NCBI(keggId, 'nucleotide', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                                 bySimilarGenes = bySimilarGenes))
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}

.getKEGG2NCBIGene <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    tryCatch(
        {
            return(.getKEGG2NCBI(keggId, 'gene', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                                 bySimilarGenes = bySimilarGenes))
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}

#########################
# KEGG database to CARD #
#########################

# Auxiliar function to retrieve intermediate NCBI ids and check their translations
# to CARD on the go
.getKEGG2CARDexhaustiveMappingAux <- function(keggId, ncbiDB, bySimilarGenes){
    aroIndex <- .loadAROIndex()

    .updateProgressBar(paste("Trying to get a direct translation for KEGG id ", keggId," to CARD", sep=""))

    # Direct translation
    translationsDT <- .getKEGG2NCBIDT(keggId)
    if(!identical(translationsDT, character(0))){
        translationsDT <- lapply(.auxNCBI2CARD(list('DT' = translationsDT),ncbiDB, aroIndex), unique)
        translationsDT <- translationsDT[lengths(translationsDT) > 0]
        if(length(translationsDT[['DT']]) > 0){
            return(translationsDT)
        }
    }

    .updateProgressBar(paste("Translating KEGG id ", keggId," to CARD through similar genes strategy", sep=""))

    # Through UniProt
    translationsTUP <- .getKEGG2NCBITUP(keggId, ncbiDB, exhaustiveMapping = TRUE, bySimilarGenes = bySimilarGenes)
    if(length(translationsTUP)>0){
        translationsTUP <- lapply(.auxNCBI2CARD(translationsTUP, ncbiDB, aroIndex), unique)
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
.getKEGG2CARD <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')

    .updateProgressBar("Searching CARD database and translating to KEGG")

    tryCatch(
        {
            .checkKEGGIdExists(keggId)

            aroIndex <- .loadAROIndex()
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
                proteinId <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getKEGG2NCBIProtein(keggId, exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes)))
                if(length(proteinId)>0){
                    result <- lapply(proteinId, FUN = function(x) unlist(sapply(x, USE.NAMES = FALSE, FUN = function(auxId){
                        auxId <- strsplit(auxId,'.', fixed = TRUE)[[1]][[1]]
                        return(aroIndex[gsub("\\..*", "", aroIndex$Protein.Accession) == auxId, "ARO.Accession"])
                    })))
                    result <- result[lengths(result) > 0]
                }

                nucleotideId <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getKEGG2NCBINucleotide(keggId, exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes)))
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
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}

########################
# Functions interfaces #
########################

getKEGG2UniProt <- function(keggId, exhaustiveMapping = FALSE) {
    .initializeProgressBar(100)
    .updateProgressBar("Translating from KEGG to UniProt")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(keggId, exhaustiveMapping = FALSE) {
                return(.retryHandler(f, keggId, exhaustiveMapping))})
        })(.getKEGG2UniProt), vectorize.args = c('keggId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(keggId, exhaustiveMapping = FALSE) {
        result <- .resultParser(keggId, exhaustiveMapping, vectorizedFunc(keggId, exhaustiveMapping))
        .terminateProgressBar()
        return(result)
    })(keggId, exhaustiveMapping))
}

getKEGG2NCBIProtein <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
    .initializeProgressBar(100)
    .updateProgressBar("Translating from KEGG to NCBI Protein")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
                return(.retryHandler(f, keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))})
        })(.getKEGG2NCBIProtein), vectorize.args = c('keggId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
        result <- .resultParser(keggId, exhaustiveMapping, vectorizedFunc(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))
        .terminateProgressBar()
        return(result)
    })(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))
}

getKEGG2NCBINucleotide <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
    .initializeProgressBar(100)
    .updateProgressBar("Translating from KEGG to NCBI Nucleotide")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
                return(.retryHandler(f, keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))})
        })(.getKEGG2NCBINucleotide), vectorize.args = c('keggId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
        result <- .resultParser(keggId, exhaustiveMapping, vectorizedFunc(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))
        .terminateProgressBar()
        return(result)
    })(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))
}

getKEGG2NCBIGene <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
    .initializeProgressBar(100)
    .updateProgressBar("Translating from KEGG to NCBI Gene")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
                return(.retryHandler(f, keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))})
        })(.getKEGG2NCBIGene), vectorize.args = c('keggId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
        result <- .resultParser(keggId, exhaustiveMapping, vectorizedFunc(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))
        .terminateProgressBar()
        return(result)
    })(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))
}

getKEGG2CARD <- function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
    .initializeProgressBar(100)
    .updateProgressBar("Translating from KEGG to CARD")

    .checkIfCARDIsDownloaded()
    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
                return(.retryHandler(f, keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))})
        })(.getKEGG2CARD), vectorize.args = c('keggId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
        result <- .resultParser(keggId, exhaustiveMapping, vectorizedFunc(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))
        .terminateProgressBar()
        return(result)
    })(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes))
}

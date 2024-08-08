###################################
# CARD database to NCBI databases #
###################################

.getCARD2NCBIAux <- function(cardId){
    .updateProgressBar(paste("Accessing CARD database and translating id ", cardId," to NCBI", sep=""))

    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- .loadAROIndex()
            resultId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]
            return(resultId)
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}

.getCARD2NCBIGene <- function(cardId, exhaustiveMapping = FALSE){
    .checkBoolean(exhaustiveMapping)
    .updateProgressBar(paste("Accessing CARD database and translating id ", cardId," to NCBI", sep=""))

    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- .loadAROIndex()
            nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
            proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession
            result <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getNCBIProtein2NCBIGene(proteinId, exhaustiveMapping = exhaustiveMapping)))

            if(exhaustiveMapping || length(result)==0){
                result <- unique(c(result, .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getNCBINucleotide2NCBIGene(nucleotideId, exhaustiveMapping = exhaustiveMapping)))))
                if(!exhaustiveMapping && length(result)>0){
                    result <- result[[1]]
                }
            }
            return(result)
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}


############################
# CARD database to UniProt #
############################

.getCARD2UniProt <- function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE){
    .checkBoolean(exhaustiveMapping)
    .checkBoolean(detailedMapping)
    .updateProgressBar(paste("Accessing CARD database and translating id ", cardId," to Uniprot", sep=""))

    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- .loadAROIndex()
            nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
            proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession

            result <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getNCBIProtein2UniProt(proteinId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = TRUE)))

            if(exhaustiveMapping || length(result)==0){
                resultAux <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getNCBINucleotide2UniProt(nucleotideId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = TRUE)))
                if(!detailedMapping){
                    result <- unique(c(result, resultAux))
                    if(is.null(result)){result <- character(0)}
                }else{
                    result <- .mergeNamedLists(result, resultAux)
                    result <- lapply(result, unique)
                    result <- result[lengths(result) > 0]
                    if(length(result)==0){result <- list()}
                }
            }

            return(result)
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}

#########################
# CARD database to KEGG #
#########################

.getCARD2KEGG <- function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(byIdenticalProteins, 'byIdenticalProteins')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    .updateProgressBar(paste("Accessing CARD database and translating id ", cardId," to KEGG", sep=""))

    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- .loadAROIndex()
            nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
            proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession

            result <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getNCBIProtein2KEGG(proteinId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes)))

            if(exhaustiveMapping || length(result)==0){
                aux <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getNCBINucleotide2KEGG(nucleotideId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes)))
                if(detailedMapping){
                    result <- lapply(.mergeNamedLists(result, aux), unique)
                }else{
                    result <- unique(c(result, aux))
                }
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

# Create the interface function
getCARD2NCBIProtein <- function(cardId) {
    .checkIfCARDIsDownloaded()
    # First, we initialize the progress bar
    .initializeProgressBar(100)
    .updateProgressBar("Translating from CARD to NCBI Protein")
    # Second, we vectorize the original function
    vectorizedFunc <- Vectorize(
        (function(f){
            # First, we wrap the memoised function with a layer that retries the request if errors are
            # raised and also handles them
            return(function(cardId) {
                return(.retryHandler(f,cardId)$Protein.Accession)
            })
        })(.getCARD2NCBIAux), vectorize.args = c('cardId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    # Third, we return the vectorised function, expliciting the original function parameters', default
    # values and parsing the result (.resultParser receives the input IDs, the exhaustiveMapping variable
    # if exists, if not, a NULL, and the vectorizedFunc result)
    return((function(cardId) {
        result <- .resultParser(cardId, NULL, vectorizedFunc(cardId))
        .terminateProgressBar()
        return(result)
    })(cardId))
}

getCARD2NCBINucleotide <- function(cardId) {
    .checkIfCARDIsDownloaded()
    .initializeProgressBar(100)
    .updateProgressBar("Translating from CARD to NCBI Nucleotide")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(cardId) {
                return(.retryHandler(f, cardId)$DNA.Accession)})
        })(.getCARD2NCBIAux), vectorize.args = c('cardId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(cardId) {
        result <- .resultParser(cardId, NULL, vectorizedFunc(cardId))
        .terminateProgressBar()
        return(result)
    })(cardId))
}

getCARD2NCBIGene <- function(cardId, exhaustiveMapping = FALSE) {
    .checkIfCARDIsDownloaded()
    .initializeProgressBar(100)
    .updateProgressBar("Translating from CARD to NCBI Gene")
    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(cardId, exhaustiveMapping = FALSE) {
                return(.retryHandler(f, cardId, exhaustiveMapping))})
        })(.getCARD2NCBIGene), vectorize.args = c('cardId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(cardId, exhaustiveMapping = FALSE) {
        result <- .resultParser(cardId, exhaustiveMapping, vectorizedFunc(cardId, exhaustiveMapping))
        .terminateProgressBar()
        return(result)
    })(cardId, exhaustiveMapping))
}

getCARD2UniProt <- function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE) {
    .checkIfCARDIsDownloaded()
    .initializeProgressBar(100)
    .updateProgressBar("Translating from CARD to UniProt")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE) {
                return(.retryHandler(f, cardId, exhaustiveMapping, detailedMapping))})
        })(.getCARD2UniProt), vectorize.args = c('cardId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)

    return((function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE) {
        result <- .resultParser(cardId, exhaustiveMapping, vectorizedFunc(cardId, exhaustiveMapping, detailedMapping))
        .terminateProgressBar()
        return(result)
    })(cardId, exhaustiveMapping, detailedMapping))
}

getCARD2KEGG <- function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE) {
    .checkIfCARDIsDownloaded()
    .initializeProgressBar(100)
    .updateProgressBar("Translating from CARD to KEGG")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE) {
                return(.retryHandler(f, cardId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes))})
        })(.getCARD2KEGG), vectorize.args = c('cardId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)

    return((function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE) {
        result <-.resultParser(cardId, exhaustiveMapping, vectorizedFunc(cardId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes))
        .terminateProgressBar()
        return(result)
    })(cardId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes))

}

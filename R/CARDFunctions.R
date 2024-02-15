source("R/utilsFunctions.R")

###################################
# CARD database to NCBI databases #
###################################

.getCARD2NCBIAux <- function(cardId){
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
    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- .loadAROIndex()
            nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
            proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession
            result <- suppressMessages(.cacheErrorHandler(raiseError=TRUE, .getNCBIProtein2NCBIGene(proteinId, exhaustiveMapping = exhaustiveMapping)))

            if(exhaustiveMapping | length(result)==0){
                result <- unique(c(result, suppressMessages(.cacheErrorHandler(raiseError=TRUE, .getNCBINucleotide2NCBIGene(nucleotideId, exhaustiveMapping = exhaustiveMapping)))))
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
    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- .loadAROIndex()
            nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
            proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession

            result <- suppressMessages(.cacheErrorHandler(raiseError=TRUE, .getNCBIProtein2UniProt(proteinId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = TRUE)))

            if(exhaustiveMapping | length(result)==0){
                resultAux <- suppressMessages(.cacheErrorHandler(raiseError=TRUE, .getNCBINucleotide2UniProt(nucleotideId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = TRUE)))
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

    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- .loadAROIndex()
            nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
            proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession

            result <- suppressMessages(.cacheErrorHandler(raiseError=TRUE, .getNCBIProtein2KEGG(proteinId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes)))

            if(exhaustiveMapping | length(result)==0){
                aux <- suppressMessages(.cacheErrorHandler(raiseError=TRUE, .getNCBINucleotide2KEGG(nucleotideId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes)))
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
        return(.resultParser(cardId, NULL, vectorizedFunc(cardId)))
    })(cardId))
}

getCARD2NCBINucleotide <- function(cardId) {
    .checkIfCARDIsDownloaded()
    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(cardId) {
                return(.retryHandler(f, cardId)$DNA.Accession)})
        })(.getCARD2NCBIAux), vectorize.args = c('cardId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(cardId) {
        return(.resultParser(cardId, NULL, vectorizedFunc(cardId)))
    })(cardId))
}

getCARD2NCBIGene <- function(cardId, exhaustiveMapping = FALSE) {
    .checkIfCARDIsDownloaded()
    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(cardId, exhaustiveMapping = FALSE) {
                return(.retryHandler(f, cardId, exhaustiveMapping))})
        })(.getCARD2NCBIGene), vectorize.args = c('cardId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(cardId, exhaustiveMapping = FALSE) {
        return(.resultParser(cardId, exhaustiveMapping, vectorizedFunc(cardId, exhaustiveMapping)))
    })(cardId, exhaustiveMapping))
}

getCARD2UniProt <- function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE) {
    .checkIfCARDIsDownloaded()
    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE) {
                return(.retryHandler(f, cardId, exhaustiveMapping, detailedMapping))})
        })(.getCARD2UniProt), vectorize.args = c('cardId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE) {
        return(.resultParser(cardId, exhaustiveMapping, vectorizedFunc(cardId, exhaustiveMapping, detailedMapping)))
    })(cardId, exhaustiveMapping, detailedMapping))
}

getCARD2KEGG <- function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE) {
    .checkIfCARDIsDownloaded()
    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE) {
                return(.retryHandler(f, cardId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes))})
        })(.getCARD2KEGG), vectorize.args = c('cardId'), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE) {
        return(.resultParser(cardId, exhaustiveMapping, vectorizedFunc(cardId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes)))
    })(cardId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes))
}

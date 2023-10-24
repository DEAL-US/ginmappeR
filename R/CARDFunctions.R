
###################################
# CARD database to NCBI databases #
###################################

getCARD2NCBIProtein <- function(cardId){
    .checkIfCARDIsDownloaded()
    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- read.csv(paste(options("cardPath")[[1]],'/card-data/aro_index.tsv',sep=''), sep='\t')
            proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession
            return(proteinId)
        },
        error = function(e) {
            warning(conditionMessage(e), "\n", call. = FALSE, noBreaks. = TRUE, immediate. = TRUE)
            return(NULL)
        }
    )

}

getCARD2NCBINucleotide <- function(cardId){
    .checkIfCARDIsDownloaded()
    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- read.csv(paste(options("cardPath")[[1]],'/card-data/aro_index.tsv',sep=''), sep='\t')
            nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
            return(nucleotideId)
        },
        error = function(e) {
            warning(conditionMessage(e), "\n", call. = FALSE, noBreaks. = TRUE, immediate. = TRUE)
            return(NULL)
        }
    )
}

getCARD2NCBIGene <- function(cardId, exhaustiveMapping = FALSE){
    .checkIfCARDIsDownloaded()
    .checkBoolean(exhaustiveMapping)
    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- read.csv(paste(options("cardPath")[[1]],'/card-data/aro_index.tsv',sep=''), sep='\t')
            nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
            proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession

            result <- unique(c(getNCBIProtein2NCBIGene(proteinId, exhaustiveMapping = exhaustiveMapping), getNCBINucleotide2NCBIGene(nucleotideId, exhaustiveMapping = exhaustiveMapping)))

            if(!exhaustiveMapping & length(result)>0){
                return(result[1])
            }
            return(result)
        },
        error = function(e) {
            warning(conditionMessage(e), "\n", call. = FALSE, noBreaks. = TRUE, immediate. = TRUE)
            return(NULL)
        }
    )
}

############################
# CARD database to UniProt #
############################

getCARD2UniProt <- function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE){
    .checkIfCARDIsDownloaded()
    .checkBoolean(exhaustiveMapping)
    .checkBoolean(detailedMapping)
    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- read.csv(paste(options("cardPath")[[1]],'/card-data/aro_index.tsv',sep=''), sep='\t')
            nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
            proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession

            result <- getNCBIProtein2UniProt(proteinId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = TRUE)
            if(exhaustiveMapping | length(result)==0){
                if(!detailedMapping){
                    result <- unique(c(result, getNCBINucleotide2UniProt(nucleotideId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = TRUE)))
                    if(is.null(result)){result <- character(0)}
                }else{
                    result <- .mergeNamedLists(result, getNCBINucleotide2UniProt(nucleotideId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = TRUE))
                    result <- lapply(result, unique)
                    result <- result[lengths(result) > 0]
                    if(length(result)==0){result <- list()}
                }
            }
            return(result)
        },
        error = function(e) {
            warning(conditionMessage(e), "\n", call. = FALSE, noBreaks. = TRUE, immediate. = TRUE)
            return(NULL)
        }
    )
}

#########################
# CARD database to KEGG #
#########################

getCARD2KEGG <- function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){

    .checkIfCARDIsDownloaded()
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(byIdenticalProteins, 'byIdenticalProteins')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')

    tryCatch(
        {
            .checkCARDIdExists(cardId)

            aroIndex <- read.csv(paste(options("cardPath")[[1]],'/card-data/aro_index.tsv',sep=''), sep='\t')
            nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
            proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession

            result <- getNCBIProtein2KEGG(proteinId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes)

            if(exhaustiveMapping | length(result)==0){
                aux <- getNCBINucleotide2KEGG(nucleotideId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes)
                if(detailedMapping){
                    result <- lapply(.mergeNamedLists(result, aux), unique)
                }else{
                    result <- unique(c(result, aux))
                }
            }
            return(result)
        },
        error = function(e) {
            warning(conditionMessage(e), "\n", call. = FALSE, noBreaks. = TRUE, immediate. = TRUE)
            return(NULL)
        }
    )
}





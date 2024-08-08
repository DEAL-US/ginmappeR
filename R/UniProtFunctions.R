############################
# UniProt database to KEGG #
############################

# Direct translation method
.getUniProt2KEGGDT <- function(upId){
    keggId <- keggConv("genes", sprintf('uniprot:%s', upId))
    if(length(keggId) == 0){
        return(character(0))
    }
    return(keggId[[1]])
}

# Method to retrieve similar genes by clusters of a % of identity
.getUniProtSimilarGenes <- function(upId, clusterIdentity = '1.0', clusterNames = FALSE){
    .checkBoolean(clusterNames, 'clusterNames')
    .checkClusterIdentity(clusterIdentity)

    # Auxiliar function
    getClusterTypeBySimilarity <- function(similarity){
        return(switch (similarity,
                       "1.0" = "uniref_cluster_100",
                       "0.9" = "uniref_cluster_90",
                       "0.5" = "uniref_cluster_50",
                       "Error"
        ))
    }

    tryCatch(
        {
            .checkUniProtIdExists(upId)

            .updateProgressBar(paste("Connecting to UniProt API and retrieving cluster of id ", upId, sep=""))

            # Get clusters list
            clusterSearchUrl <- sprintf('https://rest.uniprot.org/uniref/search?fields=id%%2Ccount&format=tsv&query=(%s%%20AND%%20identity=%s)&size=500', upId, clusterIdentity)
            clusterList <- NA
            while(identical(clusterList, NA)){
                try(clusterList <- read.csv(clusterSearchUrl, sep="\t"))
            }

            genesList <- list()
            # Iterate over clusters list to retrieve all similar genes
            for(i in 1:length(clusterList$Cluster.ID)){
                clusterName <- clusterList$Cluster.ID[i]
                .updateProgressBar(paste("Retrieving and parsing genes of cluster ", clusterName, sep=""))

                clusterSize <- clusterList$Size[i]
                similarGenesSearchUrl <- sprintf('https://rest.uniprot.org/uniprotkb/search?fields=accession&format=tsv&query=%%28%%28%s%%3A%s%%29%%20AND%%20NOT%%20%%28accession%%3A%s%%29%%29&size=500',
                                                 getClusterTypeBySimilarity(clusterIdentity) , clusterName, upId)

                # Handle UniProt's API 500 genes per request issue
                auxList <- character(0)
                request <- GET(similarGenesSearchUrl)
                auxList <- c(auxList, c(read.csv(text=content(request, as='text', encoding='UTF-8'), sep='\t')$Entry))

                while(!is.null(request$headers$link)){
                    link <- strsplit(strsplit(request$headers$link, '<')[[1]],'>')[[2]][[1]]
                    request <- GET(link)
                    auxList <- c(auxList, c(read.csv(text=content(request, as='text', encoding='UTF-8'), sep='\t')$Entry))
                }

                if(identical(logical(0), auxList)){auxList<-character(0)}
                genesList[[clusterName]] <- auxList
            }
            if(!clusterNames){
                genesList <- unname(unlist(genesList, recursive = FALSE))
            }
            return(genesList)
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}

# Method to translate UniProt to KEGG through Similar Genes Translation method, that is,
# translating genes from UniProt clusters of similar genes
.getUniProt2KEGGSGT <- function(upId, exhaustiveMapping = FALSE){
    result <- list()

    # Iterate over clusters identities
    for(clusterIdentity in c('1.0','0.9','0.5')){
        similarGenes <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getUniProtSimilarGenes(upId, clusterIdentity, TRUE)))
        # Iterate over clusters of the same identity if there were
        for(cluster in names(similarGenes)){
            cluster <- similarGenes[[cluster]]
            # Check whether the cluster has similar genes
            if(length(cluster)!=0){
                # For every similar gene, try to translate to KEGG
                for(i in 1:length(cluster)){
                    translation <- .getUniProt2KEGGDT(cluster[i])
                    # If there is a translation, return the identity of the cluster,
                    # the similar gene selected and the translation to KEGG
                    if(!identical(translation, character(0))){
                        if(!exhaustiveMapping){
                            result[[clusterIdentity]] <- c(translation)
                            return(result)
                        }else{
                            result[[clusterIdentity]] <- append(result[[clusterIdentity]], translation)
                        }
                    }
                }
            }
        }
    }
    # Return result list, if none found, would be empty
    return(lapply(result, unique))
}


# Main function
# Tries direct translation method and, if requested, by similar genes clusters
# method too
.getUniProt2KEGG <- function(upId, exhaustiveMapping = FALSE, bySimilarGenes = TRUE, detailedMapping = FALSE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    .checkBoolean(detailedMapping, 'detailedMapping')

    tryCatch(
        {
            .checkUniProtIdExists(upId)
            result <- list()
            translationsDT <- .getUniProt2KEGGDT(upId)

            .updateProgressBar(paste("Trying to get a direct translation for UniProt id ", upId," to KEGG", sep=""))

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

            .updateProgressBar(paste("Translating UniProt id ", upId," to KEGG through similar genes strategy", sep=""))

            # If required, try to translate using similar genes method
            if(bySimilarGenes){
                if(!exhaustiveMapping){
                    translationsSGT <- .getUniProt2KEGGSGT(upId, FALSE)
                    if(length(translationsSGT)>0){
                        result[[names(translationsSGT[1])]] <- translationsSGT[[names(translationsSGT[1])]][[1]]
                    }
                }else{
                    result <- append(result, .getUniProt2KEGGSGT(upId, TRUE))
                }
            }
            # To return a char vector of IDs if detailedMapping is false, simply extract
            # all IDs from result in a char vector and return it
            if(!detailedMapping){
                result <- unique(unlist(result, recursive = FALSE, use.names = FALSE))
                if(is.null(result)){result <- character(0)}
            }
            return(result)
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}

############################
# UniProt database to NCBI #
############################

# Direct translation through UniProt.ws mapper
.getUniProt2NCBIDT <- function(upId){
    query <- character(0)
    query <- c(query, .suppressOutput(mapUniProt(
        from='UniProtKB_AC-ID',
        to='EMBL-GenBank-DDBJ_CDS',
        query=c(upId),
        columns=c('accession')
    ))$To)
    query <- c(query, .suppressOutput(mapUniProt(
        from='UniProtKB_AC-ID',
        to='RefSeq_Protein',
        query=c(upId),
        columns=c('accession')
    ))$To)
    query <- c(query, .suppressOutput(mapUniProt(
        from='UniProtKB_AC-ID',
        to='RefSeq_Nucleotide',
        query=c(upId),
        columns=c('accession')
    ))$To)
    if(identical(logical(0),query) || identical(NULL,query)){
        return(character(0))
    }
    return(query)
}

# 2nd: using UP mapper, look in UP clusters of identity and try to convert to NCBI.
.getUniProt2NCBISGT <- function(upId, exhaustiveMapping = FALSE, ncbiDB){
    result <- list()

    # Iterate over clusters identities
    for(clusterIdentity in c('1.0','0.9','0.5')){
        similarGenes <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getUniProtSimilarGenes(upId, clusterIdentity, clusterNames = FALSE)))
        # Check whether there are similar genes
        if(length(similarGenes)!=0){
            # For every similar gene, try to translate to NCBI
            for(i in 1:length(similarGenes)){
                translation <- .getUniProt2NCBIDT(similarGenes[i])
                if(!identical(translation, character(0))){
                    # Filter translations to get only those of the desired ncbiDB
                    translation <- translation[which(sapply(translation, .checkIdInNCBIDataBase, ncbiDB=ncbiDB))]
                    # If there is a translation, return the identity of the cluster,
                    # the similar gene selected and the translation to NCBI
                    if(!identical(translation, character(0))){
                        if(!exhaustiveMapping){
                            result[[clusterIdentity]] <- c(translation)
                            return(result)
                        }else{
                            result[[clusterIdentity]] <- append(result[[clusterIdentity]], translation)
                        }
                    }
                }
            }
        }
    }
    # Return result list, if none found, would be empty
    return(lapply(result, unique))
}

# Warning: the combination of exhaustiveMapping=TRUE and bySimilarGenes=TRUE here takes a really
# long time to finish because of lots of translations through UniProt API that is slower than NCBI's
### Main function
.getUniProt2NCBI <- function(upId, ncbiDB = 'protein', exhaustiveMapping = FALSE, bySimilarGenes = TRUE, detailedMapping = FALSE){
    .checkUniProtIdExists(upId)

    result <- list()
    translationsDT <- .getUniProt2NCBIDT(upId)
    translationsDT <- translationsDT[which(sapply(translationsDT, .checkIdInNCBIDataBase, ncbiDB=ncbiDB))]

    .updateProgressBar(paste("Trying to get a direct translation for UniProt id ", upId," to NCBI", sep=""))

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

    .updateProgressBar(paste("Translating UniProt id ", upId," to NCBI through similar genes strategy", sep=""))

    # If required, try to translate using similar genes method
    if(bySimilarGenes){
        if(!exhaustiveMapping){
            translationsSGT <- .getUniProt2NCBISGT(upId, FALSE, ncbiDB = ncbiDB)
            if(length(translationsSGT)>0){
                result[[names(translationsSGT[1])]] <- translationsSGT[[names(translationsSGT[1])]][[1]]
            }
        }else{
            result <- append(result, .getUniProt2NCBISGT(upId, TRUE, ncbiDB = ncbiDB))
        }
    }

    if(!detailedMapping){
        result <- unique(unlist(result, recursive = FALSE, use.names = FALSE))
        if(is.null(result)){result <- character(0)}
    }

    return(result)
}

.getUniProt2NCBIProtein <- function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    .checkBoolean(detailedMapping, 'detailedMapping')

    tryCatch(
        {
            return(.getUniProt2NCBI(upId, 'protein', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                                    bySimilarGenes = bySimilarGenes))
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}

.getUniProt2NCBINucleotide <- function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    .checkBoolean(detailedMapping, 'detailedMapping')

    tryCatch(
        {
            return(.getUniProt2NCBI(upId, 'nucleotide', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                                    bySimilarGenes = bySimilarGenes))
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}

.getUniProt2NCBIGene <- function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    .checkBoolean(detailedMapping, 'detailedMapping')

    tryCatch(
        {
            return(.getUniProt2NCBI(upId, 'gene', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                                    bySimilarGenes = bySimilarGenes))
        },
        error = function(e) {return(.errorMessageHandler(e))},
        warning = function(e) {return(.warningMessageHandler(e))}
    )
}


###################
# UniProt to CARD #
###################

# Auxiliar function to retrieve intermediate NCBI ids and check their translations
# to CARD on the go
.getUniProt2CARDexhaustiveMappingAux <- function(upId, ncbiDB, bySimilarGenes){
    aroIndex <- .loadAROIndex()

    .updateProgressBar(paste("Trying to get a direct translation for Uniprot id ", upId," to CARD", sep=""))

    # Direct translation
    translationsDT <- .getUniProt2NCBIDT(upId)
    if(!identical(translationsDT, character(0))){
        translationsDT <- lapply(.auxNCBI2CARD(list('DT' = translationsDT),ncbiDB, aroIndex), unique)
        translationsDT <- translationsDT[lengths(translationsDT) > 0]
        if(length(translationsDT[['DT']]) > 0){
            return(translationsDT)
        }
    }
    .updateProgressBar(paste("Translating UniProt id ", upId," to CARD through similar genes strategy", sep=""))

    # Through UniProt
    if(bySimilarGenes){
        translationsTUP <- .getUniProt2NCBISGT(upId, exhaustiveMapping = TRUE,  ncbiDB)
        if(length(translationsTUP)>0){
            translationsTUP <- lapply(.auxNCBI2CARD(translationsTUP, ncbiDB, aroIndex), unique)
            translationsTUP <- translationsTUP[lengths(translationsTUP) > 0]
            if(length(translationsTUP)>0){
                return(translationsTUP)
            }else{
                return(list())
            }
        }
    }else{
        return(list())
    }

}

.getUniProt2CARD <- function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')

    .updateProgressBar(paste("Translating Uniprot id ", upId," to CARD", sep=""))

    tryCatch(
        {
            .checkUniProtIdExists(upId)

            aroIndex <- .loadAROIndex()
            result <- list()

            # Handle not exhaustiveMapping case, trying to get the fastest result present in CARD
            if(!exhaustiveMapping){
                proteinId <-.getUniProt2CARDexhaustiveMappingAux(upId, 'protein', bySimilarGenes)
                if(length(proteinId)>0){
                    proteinId <- Filter(Negate(is.null), proteinId[c("DT", "1.0", "0.9", "0.5")])
                    result[names(proteinId)[1]] <- proteinId[[1]][[1]]
                }else{
                    nucleotideId <-.getUniProt2CARDexhaustiveMappingAux(upId, 'nucleotide', bySimilarGenes)
                    if(length(nucleotideId)>0){
                        nucleotideId <- Filter(Negate(is.null), nucleotideId[c("DT", "1.0", "0.9", "0.5")])
                        result[names(nucleotideId)[1]] <- nucleotideId[[1]][[1]]
                    }
                }
            }else{
                # Handle exhaustiveMapping case, trying to retrieve all possible cases
                proteinId <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getUniProt2NCBIProtein(upId, exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes)))
                if(length(proteinId)>0){
                    result <- lapply(proteinId, FUN = function(x) unlist(sapply(x, USE.NAMES = FALSE, FUN = function(auxId){
                        auxId <- strsplit(auxId,'.', fixed = TRUE)[[1]][[1]]
                        return(aroIndex[gsub("\\..*", "", aroIndex$Protein.Accession) == auxId, "ARO.Accession"])
                    })))
                    result <- result[lengths(result) > 0]
                }

                nucleotideId <- .suppressOutput(.cacheErrorHandler(raiseError=TRUE, .getUniProt2NCBINucleotide(upId, exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes)))
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

getUniProtSimilarGenes <- function(upId, clusterIdentity = '1.0', clusterNames = FALSE) {
    .initializeProgressBar(100)
    .updateProgressBar("Accessing UniProt similar genes database")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(upId, clusterIdentity = '1.0', clusterNames = FALSE) {
                return(.retryHandler(f, upId, clusterIdentity, clusterNames))})
        })(.getUniProtSimilarGenes), vectorize.args = c("upId", "clusterIdentity"), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(upId, clusterIdentity = '1.0', clusterNames = FALSE) {
        result <- .resultParser(upId, TRUE, vectorizedFunc(upId, clusterIdentity, clusterNames))
        .terminateProgressBar()
        return(result)
    })(upId, clusterIdentity, clusterNames))
}

getUniProt2KEGG <- function(upId, exhaustiveMapping = FALSE, bySimilarGenes = TRUE, detailedMapping = FALSE) {
    .initializeProgressBar(100)
    .updateProgressBar("Translating from UniProt to KEGG")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(upId, exhaustiveMapping = FALSE, bySimilarGenes = TRUE, detailedMapping = FALSE) {
                return(.retryHandler(f, upId, exhaustiveMapping, bySimilarGenes, detailedMapping))})
        })(.getUniProt2KEGG), vectorize.args = c("upId"), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(upId, exhaustiveMapping = FALSE, bySimilarGenes = TRUE, detailedMapping = FALSE) {
        result <- .resultParser(upId, exhaustiveMapping, vectorizedFunc(upId, exhaustiveMapping, bySimilarGenes, detailedMapping))
        .terminateProgressBar()
        return(result)
    })(upId, exhaustiveMapping, bySimilarGenes, detailedMapping))
}

getUniProt2NCBIProtein <- function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
    .initializeProgressBar(100)
    .updateProgressBar("Translating from UniProt to NCBI Protein")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
                return(.retryHandler(f, upId, exhaustiveMapping, detailedMapping, bySimilarGenes))})
        })(.getUniProt2NCBIProtein), vectorize.args = c("upId"), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
        result <- .resultParser(upId, exhaustiveMapping, vectorizedFunc(upId, exhaustiveMapping, detailedMapping, bySimilarGenes))
        .terminateProgressBar()
        return(result)
    })(upId, exhaustiveMapping, detailedMapping, bySimilarGenes))
}

getUniProt2NCBINucleotide <- function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
    .initializeProgressBar(100)
    .updateProgressBar("Translating from UniProt to NCBI Nucleotide")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
                return(.retryHandler(f, upId, exhaustiveMapping, detailedMapping, bySimilarGenes))})
        })(.getUniProt2NCBINucleotide), vectorize.args = c("upId"), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
        result <- .resultParser(upId, exhaustiveMapping, vectorizedFunc(upId, exhaustiveMapping, detailedMapping, bySimilarGenes))
        .terminateProgressBar()
        return(result)
    })(upId, exhaustiveMapping, detailedMapping, bySimilarGenes))
}

getUniProt2NCBIGene <- function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
    .initializeProgressBar(100)
    .updateProgressBar("Translating from UniProt to NCBI Gene")

    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
                return(.retryHandler(f, upId, exhaustiveMapping, detailedMapping, bySimilarGenes))})
        })(.getUniProt2NCBIGene), vectorize.args = c("upId"), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
        result <- .resultParser(upId, exhaustiveMapping, vectorizedFunc(upId, exhaustiveMapping, detailedMapping, bySimilarGenes))
        .terminateProgressBar()
        return(result)
    })(upId, exhaustiveMapping, detailedMapping, bySimilarGenes))
}

getUniProt2CARD <- function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
    .initializeProgressBar(100)
    .updateProgressBar("Translating from UniProt to CARD")

    .checkIfCARDIsDownloaded()
    vectorizedFunc <- Vectorize(
        (function(f){
            return(function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
                return(.retryHandler(f, upId, exhaustiveMapping, detailedMapping, bySimilarGenes))})
        })(.getUniProt2CARD), vectorize.args = c("upId"), USE.NAMES = FALSE, SIMPLIFY = FALSE)
    return((function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE) {
        result <- .resultParser(upId, exhaustiveMapping, vectorizedFunc(upId, exhaustiveMapping, detailedMapping, bySimilarGenes))
        .terminateProgressBar()
        return(result)
    })(upId, exhaustiveMapping, detailedMapping, bySimilarGenes))
}


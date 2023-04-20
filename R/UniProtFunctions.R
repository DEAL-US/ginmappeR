############################
# UniProt database to KEGG #
############################

# Direct translation method
.getUniProt2KEGGDT <- function(upId){
    keggId <- paste(keggConv("genes", sprintf('uniprot:%s', upId)), collapse = ';')
    if(identical(keggId, character(0))|identical(keggId, "")){
        return(character(0))
    }else{
        return(c(strsplit(keggId,';',fixed = TRUE)[[1]]))
    }
}

# Method to retrieve similar genes by clusters of a % of identity
getUniProtSimilarGenes <- function(upId, clusterIdentity = '1.0', clusterNames = FALSE){
    .checkBoolean(clusterNames, 'clusterNames')

    # Auxiliar function
    getClusterTypeBySimilarity <- function(similarity){
        return(switch (similarity,
                       "1.0" = "uniref_cluster_100",
                       "0.9" = "uniref_cluster_90",
                       "0.5" = "uniref_cluster_50",
                       "Error"
        ))
    }

    .checkClusterIdentity(clusterIdentity)
    .checkUniProtIdExists(upId)

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
}

# Method to translate UniProt to KEGG through Similar Genes Translation method, that is,
# translating genes from UniProt clusters of similar genes
.getUniProt2KEGGSGT <- function(upId, exhaustiveMapping = FALSE){
    result <- list()

    # Iterate over clusters identities
    for(clusterIdentity in c('1.0','0.9','0.5')){
        similarGenes <- getUniProtSimilarGenes(upId, clusterIdentity, TRUE)
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
getUniProt2KEGG <- function(upId, exhaustiveMapping = FALSE, bySimilarGenes = TRUE, detailedMapping = FALSE){
    .checkUniProtIdExists(upId)
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    .checkBoolean(detailedMapping, 'detailedMapping')

    result <- list()
    translationsDT <- .getUniProt2KEGGDT(upId)

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
}

############################
# UniProt database to NCBI #
############################

# Direct translation through UniProt.ws mapper
.getUniProt2NCBIDT <- function(upId){
    query <- character(0)
    query <- c(query, suppressWarnings(mapUniProt(
        from='UniProtKB_AC-ID',
        to='EMBL-GenBank-DDBJ_CDS',
        query=c(upId),
        columns=c('accession')
    ))$To)
    query <- c(query, suppressWarnings(mapUniProt(
        from='UniProtKB_AC-ID',
        to='RefSeq_Protein',
        query=c(upId),
        columns=c('accession')
    ))$To)
    query <- c(query, suppressWarnings(mapUniProt(
        from='UniProtKB_AC-ID',
        to='RefSeq_Nucleotide',
        query=c(upId),
        columns=c('accession')
    ))$To)
    if(identical(logical(0),query)|identical(NULL,query)){
        return(character(0))
    }
    return(query)
}

# 2nd: using UP mapper, look in UP clusters of identity and try to convert to NCBI.
.getUniProt2NCBISGT <- function(upId, exhaustiveMapping = FALSE, ncbiDB){
    result <- list()

    # Iterate over clusters identities
    for(clusterIdentity in c('1.0','0.9','0.5')){
        similarGenes <- getUniProtSimilarGenes(upId, clusterIdentity, clusterNames = FALSE)
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
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')
    .checkBoolean(detailedMapping, 'detailedMapping')

    result <- list()
    translationsDT <- .getUniProt2NCBIDT(upId)
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

getUniProt2NCBIProtein <- function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    return(.getUniProt2NCBI(upId, 'protein', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                         bySimilarGenes = bySimilarGenes))
}

getUniProt2NCBINucleotide <- function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    return(.getUniProt2NCBI(upId, 'nucleotide', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                            bySimilarGenes = bySimilarGenes))
}

getUniProt2NCBIGene <- function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    return(.getUniProt2NCBI(upId, 'gene', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                            bySimilarGenes = bySimilarGenes))
}



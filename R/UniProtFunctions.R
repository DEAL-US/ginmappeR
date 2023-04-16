############################
# UniProt database to KEGG #
############################

# Direct translation method
.getUniProt2KEGGDT <- function(upId){
    keggId <- paste(keggConv("genes", sprintf('uniprot:%s', upId)), collapse = ';')
    if(identical(keggId, character(0))|identical(keggId, "")){
        return(list())
    }else{
        return(as.list(strsplit(keggId,';')[[1]]))
    }
}

# Method to retrieve similar genes by clusters of a % of identity
getUniProtSimilarGenes <- function(upId, clusterIdentity = '1.0'){

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
        auxList <- list()
        request <- GET(similarGenesSearchUrl)
        auxList <- append(auxList, as.list(read.csv(text=content(request, as='text', encoding='UTF-8'), sep='\t')$Entry))

        while(!is.null(request$headers$link)){
            link <- strsplit(strsplit(request$headers$link, '<')[[1]],'>')[[2]][[1]]
            request <- GET(link)
            auxList <- append(auxList, as.list(read.csv(text=content(request, as='text', encoding='UTF-8'), sep='\t')$Entry))
        }

        genesList[[clusterName]] <- auxList
        return(genesList)
    }
}

# Method to translate UniProt to KEGG through Similar Genes Translation method, that is,
# translating genes from UniProt clusters of similar genes
.getUniProt2KEGGSGT <- function(upId){

    # Iterate over clusters identities
    for(clusterIdentity in c('1.0','0.9','0.5')){
        similarGenes <- getUniProtSimilarGenes(upId, clusterIdentity)
        # Iterate over clusters of the same identity if there were
        for(cluster in names(similarGenes)){
            cluster <- similarGenes[[cluster]]
            # Check whether the cluster has similar genes
            if(length(cluster)!=0){
                # For every similar gene, try to translate to KEGG
                for(i in 1:length(cluster)){
                    translation <- .getUniProt2KEGGDT(cluster[[i]])
                    # If there is a translation, return the identity of the cluster,
                    # the similar gene selected and the translation to KEGG
                    if(!identical(translation, list())){
                        return(c(clusterIdentity, cluster[[i]], translation))
                    }
                }
            }
        }
    }
    # If a translation has not been found, return an empty list
    return(list())
}

# Main function --- Explanations below in order to facilitate documentation writing
# Tries direct translation method and, if requested, by similar genes clusters
# method too
# If SGT method is selected, returns a list with two lists, one with de DT method
# result and another one with the SGT method result. For example:
# list( 'DT' = list(translation), 'SGT' = list(clusterIdentity, intermediateUPID, translation))
# list( 'DT' = list(ag:AEX08599), 'SGT' = list(0.9, C7C422, ag:CAZ39946))
getUniProt2KEGG <- function(upId, bySimilarGenes = TRUE){
    .checkUniProtIdExists(upId)

    translations <- list('DT'=.getUniProt2KEGGDT(upId))

    if(bySimilarGenes){
        translations <- append(translations, list('SGT'=.getUniProt2KEGGSGT(upId)))
    }

    return(translations)
}






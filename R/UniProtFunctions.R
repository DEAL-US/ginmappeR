############################
# UniProt database to KEGG #
############################

# Main function
getUniProt2KEGG <- function(upId){
    .checkUniProtIdExists(upId)

    #TODO
}

# Direct translation method
.getUniProt2KEGGDT <- function(upId){
    keggId <- paste(keggConv("genes", sprintf('uniprot:%s', upId)), collapse = ';')
    if(identical(keggId, character(0))|identical(keggId, "")){
        return(list())
    }else{
        return(as.list(strsplit(keggId,';')[[1]]))
    }
}

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


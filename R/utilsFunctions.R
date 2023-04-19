.checkNotNull <- function(variable, message="The given ID is "){
    if(identical(variable, "")){
        stop(paste(message, "empty", sep=""))
    }else if(is.null(variable)){
        stop(paste(message, "NULL", sep=""))
    }else if(is.nan(variable)){
        stop(paste(message, "NaN", sep=""))
    }else if(is.na(variable)){
        stop(paste(message, "NA", sep=""))
    }
}

.checkNCBIProteinIdExists <- function(proteinId){
    .checkNotNull(proteinId, 'The given NCBI Protein ID is ')
    # Example id: ADH03009.1
    # [ACCN] means we look for an Accesion Number OR [UID] Unique number assigned to publication
    invisible(capture.output(query <- entrez_search(db='protein', term=paste(proteinId,'[ACCN] OR ', proteinId, '[UID]',sep=''))))
    if(!length(query[['ids']])>0) {
        stop(paste('The given ID "', proteinId, '" is not registered in NCBI Protein database', sep = ''))
    }
}

.checkNCBINucleotideIdExists <- function(nucleotideId){
    .checkNotNull(nucleotideId, 'The given NCBI Nucleotide ID is ')
    # Example id: HM036080.1
    invisible(capture.output(query <- entrez_search(db='nucleotide', term=paste(nucleotideId,'[ACCN] OR ', nucleotideId, '[UID]',sep=''))))
    if(!length(query[['ids']])>0) {
        stop(paste('The given ID "', nucleotideId, '" is not registered in NCBI Nucleotide database', sep = ''))
    }
}

.checkNCBIGeneIdExists <- function(geneId){
    .checkNotNull(geneId, 'The given NCBI Gene ID is ')
    # Example id: WP_001082319
    invisible(capture.output(query <- entrez_search(db='gene', term=paste(geneId,'[ACCN] OR ', geneId, '[UID]',sep=''))))
    if(!length(query[['ids']])>0) {
        stop(paste('The given ID "', geneId, '" is not registered in NCBI Gene database', sep = ''))
    }
}

.checkUniProtIdExists <- function(upId){
    .checkNotNull(upId, 'The given UniProt Gene ID is ')
    # Example id not registered: P0ZTH2
    # Example id not valid: P0ZTUH2
    # Example valid id: P0DTH6
    output <- invisible(capture.output(queryUniProt(
        query = paste('accession:', upId, sep=''),
        fields = c('accession'),
        collapse = ' OR '
    )))
    if(grepl('invalid format', output[2], fixed = TRUE)){
        stop(paste('The given ID "', upId, '" has not a valid UniProt database accesion format', sep = ''))
    }else if(grepl('0 rows', output[2], fixed = TRUE)){
        stop(paste('The given ID "', upId, '" is not registered in UniProt database', sep = ''))
    }
}

.checkClusterIdentity <- function(clusterIdentity){
    if(!clusterIdentity %in% c('1.0','0.9','0.5')){
        stop(paste('The requested cluster identity "', clusterIdentity, '" is not valid', sep = ''))
    }
}

.checkBoolean <- function(value, name){
    if(!value %in% c(TRUE, FALSE)){
        stop(paste(name, ' variable value must be TRUE or FALSE, value given: ', value, sep = ''))
    }
}

.mergeNamedLists <- function(list1,list2){
    for (key in names(list2)) {
        if (key %in% names(list1)) {
            list1[[key]] <- c(list1[[key]], list2[[key]])
        } else {
            list1[[key]] <- list2[[key]]
        }
    }
    return(list1)
}

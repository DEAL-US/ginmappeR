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
    .checkNotNull(upId, 'The given UniProt ID is ')
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

.checkKEGGIdExists <- function(keggId){
    .checkNotNull(keggId, 'The given KEGG ID is ')

    # Search ID in KEGG genes database
    query <- keggFind('genes',keggId)

    # Check if there is a coincidence and if it corresponds to an unique gene
    if(!length(query)==1){
        stop(paste('The given ID "', keggId, '" is not registered in KEGG genes database', sep = ''))
    }
}

.checkCARDIdExists <- function(cardId){
    .checkNotNull(cardId, 'The given CARD ARO accession is ')

    checkedCardId <- if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)

    aroIndex <- read.csv(paste(getwd(),'/card-data/aro_index.tsv',sep=''), sep='\t')

    if(!checkedCardId %in% aroIndex$ARO.Accession) {
        stop(paste('The given ARO accession "', cardId, '" is not registered in CARD database', sep = ''))
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

.checkIdInNCBIDataBase <- function(ncbiId, ncbiDB) {
    result<- tryCatch(
        switch (ncbiDB,
                'protein' = .checkNCBIProteinIdExists(ncbiId),
                'nucleotide' = .checkNCBINucleotideIdExists(ncbiId),
                'gene' = .checkNCBIGeneIdExists(ncbiId),
        ),
        error = function(e){return(FALSE)}
    )
    if(is.null(result)){result <- TRUE}
    return(result)
}

#####################################
# CARD database retrieval functions #
#####################################

.downloadAndExtractCARD <- function(){

    pkgInfo <- get_pkg_info("ginmappeR")
    localRelativeFilenames <- c(paste('card-data.tar.bz2-', format(Sys.time(), "%a %m-%d-%Y %H-%M-%S"),sep=''))
    url <- c('https://card.mcmaster.ca/latest/data')
    # Download CARD zip if not present
    ensure_files_available(pkgInfo, localRelativeFilenames, url)

    # If it has not been extracted, extract it
    if(!file.exists(paste(getwd(),'/card-data/aro_index.tsv',sep=''))){
        zipF<- paste(get_cache_dir(pkgInfo),'/',localRelativeFilenames[1], sep='')

        # Unzip the compressed file
        untar(zipF, exdir = paste(getwd(),'/card-data',sep=''))
    }
}

updateCARDDataBase <- function(){
    message('Updating CARD database data...')
    erase_file_cache(get_pkg_info("ginmappeR"))
    unlink(paste(getwd(),'/card-data',sep=''), recursive = TRUE)
    .downloadAndExtractCARD()
    message(paste('CARD database downloaded successfully!\nLocated at ', paste(getwd(),'/card-data',sep=''), sep=''))
}

.checkIfCARDIsDownloaded <- function(){
    pkgInfo <- get_pkg_info("ginmappeR")
    if(!file.exists(paste(getwd(),'/card-data/aro_index.tsv',sep=''))){
        updateCARDDataBase()
    }else{
        message(sprintf('Using a CARD database version downloaded on %s, please consider updating it with updateCARDDataBase() function.',
                        strsplit(list.files(get_cache_dir(get_pkg_info('ginmappeR')))[[1]],'bz2-')[[1]][[2]]))
    }
}

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
    invisible(capture.output(query <- entrez_search(db='protein', term=proteinId)))
    if(!length(query[['ids']])>0) {
        stop(paste('The given ID "', proteinId, '" is not registered in NCBI Protein database', sep = ""))
    }
}

.checkNCBINucleotideIdExists <- function(nucleotideId){
    .checkNotNull(nucleotideId, 'The given NCBI Nucleotide ID is ')
    # Example id: HM036080.1
    invisible(capture.output(query <- entrez_search(db='nucleotide', term=nucleotideId)))
    if(!length(query[['ids']])>0) {
        stop(paste('The given ID "', nucleotideId, '" is not registered in NCBI Nucleotide database', sep = ""))
    }
}

.checkNCBIGeneIdExists <- function(geneId){
    .checkNotNull(geneId, 'The given NCBI Gene ID is ')
    # Example id: WP_001082319
    invisible(capture.output(query <- entrez_search(db='gene', term=geneId)))
    if(!length(query[['ids']])>0) {
        stop(paste('The given ID "', geneId, '" is not registered in NCBI Gene database', sep = ""))
    }
}




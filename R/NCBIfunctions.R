#####################################
# NCBI databases inter-translations #
#####################################

.getNCBIDatabasesLinks <- function(dbFrom="protein", id, dbTo="gene"){
    return(entrez_link(dbfrom=dbFrom, id=id, db=dbTo))
}

getNCBIGene2NCBIProtein <- function(id){
    .checkNCBIGeneIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='gene', id=id, dbTo='protein')
    result <- list()
    # Handle multiple protein IDs case
    if(!identical(query[['links']][['gene_protein']], NULL)){
        for(query_id in query[['links']][['gene_protein']]){
            protein_xml <- entrez_fetch(db = "protein", id = query_id, rettype = "xml")
            result <- append(result, xpathSApply(xmlParse(protein_xml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
        }
    }
    # Handle no possible translations case
    if(length(result)==0){
        stop("The given NCBI Gene ID has no translation to NCBI Protein database")
    }
    return(result)
}

getNCBIProtein2NCBIGene <- function(id){
    .checkNCBIProteinIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='protein', id=id, dbTo='gene')
    # Handle no possible translations case
    if(length(query[['links']][['protein_gene']])==0){
        stop("The given NCBI Protein ID has no translation to NCBI Gene database")
    }
    return(as.list(query[['links']][['protein_gene']]))
}

getNCBIProtein2NCBINucleotide <- function(id){
    .checkNCBIProteinIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='protein', id=id, dbTo='nucleotide')
    result <- list()
    # Handle multiple protein IDs case
    if(!identical(query[['links']][['protein_nuccore']], NULL)){
        for(query_id in query[['links']][['protein_nuccore']]){
            protein_xml <- entrez_fetch(db = "nucleotide", id = query_id, rettype = "xml")
            result <- append(result, xpathSApply(xmlParse(protein_xml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
        }
    }
    # Handle no possible translations case
    if(length(result)==0){
        stop("The given NCBI Protein ID has no translation to NCBI Nucleotide database")
    }
    return(result)

}

getNCBINucleotide2NCBIProtein <- function(id){
    .checkNCBINucleotideIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='nucleotide', id=id, dbTo='protein')
    result <- list()
    # Handle multiple protein IDs case
    if(!identical(query[['links']][['nuccore_protein']], NULL)){
        for(query_id in query[['links']][['nuccore_protein']]){
            protein_xml <- entrez_fetch(db = "protein", id = query_id, rettype = "xml")
            result <- append(result, xpathSApply(xmlParse(protein_xml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
        }
    }
    # Handle no possible translations case
    if(length(result)==0){
        stop("The given NCBI Nucleotide ID has no translation to NCBI Protein database")
    }
    return(result)
}

getNCBIGene2NCBINucleotide <- function(id){
    .checkNCBIGeneIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='gene', id=id, dbTo='nucleotide')
    result <- list()
    # Handle multiple protein IDs case
    if(!identical(query[['links']][['gene_nuccore']], NULL)){
        for(query_id in query[['links']][['gene_nuccore']]){
            protein_xml <- entrez_fetch(db = "nucleotide", id = query_id, rettype = "xml")
            result <- append(result, xpathSApply(xmlParse(protein_xml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
        }
    }
    # Handle no possible translations case
    if(length(result)==0){
        stop("The given NCBI Gene ID has no translation to NCBI Nucleotide database")
    }
    return(result)
}

getNCBINucleotide2NCBIGene <- function(id){
    .checkNCBINucleotideIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='nucleotide', id=id, dbTo='gene')
    # Handle no possible translations case
    if(length(query[['links']][['nuccore_gene']])==0){
        stop("The given NCBI Nucleotide ID has no translation to NCBI Gene database")
    }
    return(as.list(query[['links']][['nuccore_gene']]))
}

#####################################
# NCBI databases to UniProt #
#####################################

# Direct translation method
.getNCBI2UniProtDT <- function(ncbiId, database='EMBL-GenBank-DDBJ_CDS'){
    query <- NULL
    tryCatch(
        {
        invisible(capture.output(query <- mapUniProt(
            from=database,
            to='UniProtKB',
            query=c(ncbiId),
            columns=c('accession')
            )))
        },
        warning=function(w){}
    )
    return(query$Entry)
}

# Identical NCBI proteins
getNCBIIdenticalProteins <- function(ncbiId, format = 'ids'){
    if(!format %in% c('ids','dataframe')){
        stop('Incorrect format requested')
    }

    # Retrieve intermediate ID
    auxId <- NA
    while(identical(auxId,NA)){
        auxId <- entrez_search(db = "ipg", term=paste(ncbiId,'[ACCN] OR ', ncbiId, '[UID]',sep=''))
    }

    # Check if intermediate ID has been retrieved (shows us if there are or not identical proteins registered)
    if(length(auxId$ids)==0){
        stop('The given NCBI ID has no similar proteins registered')
    }else{
        auxId <- auxId$ids[[1]]
        identicalProteins <- NA
        while(identical(identicalProteins,NA)){
            identicalProteins <- entrez_fetch(db = 'ipg', id = auxId, rettype = 'native')
        }
        identicalProteins <- read.csv(text=identicalProteins, sep = '\t')
        identicalProteins <- identicalProteins[which(identicalProteins$Source=='RefSeq' | identicalProteins$Source=='INSDC'),]
        switch(format,
               'ids'= return(as.list(identicalProteins$Protein)),
               'dataframe' = return(identicalProteins))
    }
}


#
# # Identical proteins translation method
# .getNCBI2UniProtIP <- function(ncbiId){
#     return(ncbiId)
# }
#
# # NCBI GenBank to UniProt translation
# getNCBIGenBank2UniProt <- function(ncbiId){
#     .checkNotNull(ncbiId, 'The NCBI id given is NULL')
#     .checkNCBIIdExists(ncbiId)
#
#     # First translation strategy
#     translatedId <- .getNCBI2UniProtDT(ncbiId)
#
#     # Second translation strategy
#     if(identical(logical(0), translatedId)){
#         translatedId <- .getNCBI2UniProtIP(ncbiId)
#     }
#     return(translatedId)
# }

# TODO: en los métodos de NCBIxxx a UniProt, primero hay que comprobar que el ID exista en NCBIxxx



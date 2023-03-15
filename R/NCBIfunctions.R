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
        warning("The given NCBI Gene ID has no translation to NCBI Protein database")
    }
    return(result)
}

getNCBIProtein2NCBIGene <- function(id){
    .checkNCBIProteinIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='protein', id=id, dbTo='gene')
    # Handle no possible translations case
    if(length(query[['links']][['protein_gene']])==0){
        warning("The given NCBI Protein ID has no translation to NCBI Gene database")
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
        warning("The given NCBI Protein ID has no translation to NCBI Nucleotide database")
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
        warning("The given NCBI Nucleotide ID has no translation to NCBI Protein database")
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
        warning("The given NCBI Gene ID has no translation to NCBI Nucleotide database")
    }
    return(result)
}

getNCBINucleotide2NCBIGene <- function(id){
    .checkNCBINucleotideIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='nucleotide', id=id, dbTo='gene')
    # Handle no possible translations case
    if(length(query[['links']][['nuccore_gene']])==0){
        warning("The given NCBI Nucleotide ID has no translation to NCBI Gene database")
    }
    return(as.list(query[['links']][['nuccore_gene']]))
}

#####################################
# NCBI databases to UniProt #
#####################################

# # Direct translation method
# .getNCBI2UniProtDT <- function(ncbiId, database='EMBL-GenBank-DDBJ_CDS'){
#
#     return(suppressWarnings(mapUniProt(
#         from=database,
#         to='UniProtKB',
#         query=c(ncbiId),
#         columns=c('accession')
#     ))$Entry)
# }
#
# # Identical NCBI proteins
# getNCBIIdenticalProteins <- function(ncbiId){
#     return(NULL)
# }
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


############################
# KEGG database to UniProt #
############################

getKEGG2UniProt <- function(keggId){
    .checkKEGGIdExists(keggId)

    query <- keggConv('uniprot', keggId)
    if(identical(query, character(0))|identical(query, "")){
        return(character(0))
    }else{
        return(c(strsplit(query,':',fixed = TRUE)[[1]][[2]]))
    }
}

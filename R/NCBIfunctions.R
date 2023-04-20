#####################################
# NCBI databases inter-translations #
#####################################

.getNCBIDatabasesLinks <- function(dbFrom="protein", id, dbTo="gene"){
    return(entrez_link(dbfrom=dbFrom, id=id, db=dbTo))
}

getNCBIGene2NCBIProtein <- function(id, exhaustiveMapping = FALSE){
    .checkNCBIGeneIdExists(id)
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')

    query <- .getNCBIDatabasesLinks(dbFrom='gene', id=id, dbTo='protein')
    result <- character(0)
    # Handle multiple protein IDs case
    if(!identical(query[['links']][['gene_protein']], NULL)){
        for(queryId in query[['links']][['gene_protein']]){
            proteinXml <- entrez_fetch(db = "protein", id = queryId, rettype = "xml")
            result <- c(result, xpathSApply(xmlParse(proteinXml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
            if(!exhaustiveMapping){
                return(result)
            }
        }
    }
    return(result)
}

getNCBIProtein2NCBIGene <- function(id, exhaustiveMapping = FALSE){
    .checkNCBIProteinIdExists(id)
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')

    query <- .getNCBIDatabasesLinks(dbFrom='protein', id=id, dbTo='gene')
    if(!identical(query[['links']][['protein_gene']], NULL)){
        if(!exhaustiveMapping){
            return(c(query[['links']][['protein_gene']])[1])
        }
        return(c(query[['links']][['protein_gene']]))
    }else{
        return(character(0))
    }
}

getNCBIProtein2NCBINucleotide <- function(id, exhaustiveMapping = FALSE){
    .checkNCBIProteinIdExists(id)
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')

    query <- .getNCBIDatabasesLinks(dbFrom='protein', id=id, dbTo='nucleotide')
    result <- character(0)
    # Handle multiple protein IDs case
    if(!identical(query[['links']][['protein_nuccore']], NULL)){
        for(queryId in query[['links']][['protein_nuccore']]){
            proteinXml <- entrez_fetch(db = "nucleotide", id = queryId, rettype = "xml")
            result <- c(result, xpathSApply(xmlParse(proteinXml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
            if(!exhaustiveMapping){
                return(result)
            }
        }
    }
    return(result)
}

getNCBINucleotide2NCBIProtein <- function(id, exhaustiveMapping = FALSE){
    .checkNCBINucleotideIdExists(id)
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')

    query <- .getNCBIDatabasesLinks(dbFrom='nucleotide', id=id, dbTo='protein')
    result <- character(0)
    # Handle multiple protein IDs case
    if(!identical(query[['links']][['nuccore_protein']], NULL)){
        for(queryId in query[['links']][['nuccore_protein']]){
            proteinXml <- entrez_fetch(db = "protein", id = queryId, rettype = "xml")
            result <- c(result, xpathSApply(xmlParse(proteinXml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
            if(!exhaustiveMapping){
                return(result)
            }
        }
    }
    return(result)
}

getNCBIGene2NCBINucleotide <- function(id, exhaustiveMapping = FALSE){
    .checkNCBIGeneIdExists(id)
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')

    query <- .getNCBIDatabasesLinks(dbFrom='gene', id=id, dbTo='nucleotide')
    result <- character(0)
    # Handle multiple protein IDs case
    if(!identical(query[['links']][['gene_nuccore']], NULL)){
        for(queryId in query[['links']][['gene_nuccore']]){
            proteinXml <- entrez_fetch(db = "nucleotide", id = queryId, rettype = "xml")
            result <- c(result, xpathSApply(xmlParse(proteinXml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
            if(!exhaustiveMapping){
                return(result)
            }
        }
    }
    return(result)
}

getNCBINucleotide2NCBIGene <- function(id, exhaustiveMapping = FALSE){
    .checkNCBINucleotideIdExists(id)
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')

    query <- .getNCBIDatabasesLinks(dbFrom='nucleotide', id=id, dbTo='gene')

    if(!identical(query[['links']][['nuccore_gene']], NULL)){
        if(!exhaustiveMapping){
            return(c(query[['links']][['nuccore_gene']])[1])
        }
        return(c(query[['links']][['nuccore_gene']]))
    }else{
        return(character(0))
    }
}

#############################
# NCBI databases to UniProt #
#############################

# Direct translation method
.getNCBI2UniProtDT <- function(ncbiId){
    query <- character(0)
    query <- c(query, suppressWarnings(mapUniProt(
        from='EMBL-GenBank-DDBJ_CDS',
        to='UniProtKB',
        query=c(ncbiId),
        columns=c('accession')
    ))$Entry)
    query <- c(query, suppressWarnings(mapUniProt(
        from='RefSeq_Protein',
        to='UniProtKB',
        query=c(ncbiId),
        columns=c('accession')
    ))$Entry)
    query <- c(query, suppressWarnings(mapUniProt(
        from='RefSeq_Nucleotide',
        to='UniProtKB',
        query=c(ncbiId),
        columns=c('accession')
    ))$Entry)
    if(identical(logical(0),query)|identical(NULL,query)){
        return(character(0))
    }
    return(query)
}

# Get identical NCBI proteins
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
        switch(format,
               'ids'= return(character(0)),
               'dataframe' = return(data.frame()))
        # stop('The given NCBI ID has no similar proteins registered')
    }else{
        auxId <- auxId$ids[[1]]
        identicalProteins <- NA
        while(identical(identicalProteins,NA)){
            identicalProteins <- entrez_fetch(db = 'ipg', id = auxId, rettype = 'native')
        }
        identicalProteins <- read.csv(text=identicalProteins, sep = '\t')
        identicalProteins <- identicalProteins[which(identicalProteins$Source=='RefSeq' | identicalProteins$Source=='INSDC'),]
        switch(format,
               'ids'= return(c(identicalProteins$Protein)),
               'dataframe' = return(identicalProteins))
    }
}

# Translate a batch of NCBI IDs to UniProt
.getNCBI2UniProtBatch <- function(proteins, database){

    makeQuery <- function(proteinsQuery, database){
        query <- suppressWarnings(mapUniProt(
            from=database,
            to='UniProtKB',
            query=c(proteinsQuery),
            columns=c('accession')
        ))
        return(query)
    }

    resultDf <- NA

    chunk <- 99999
    n <- nrow(proteins)
    r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
    d <- split(proteins,r)

    for(i in 1:length(d)){
        if(i==1){
            resultDf <- makeQuery(unique(d[[i]]$Protein), database)
        }else{
            resultDf <- rbind(resultDf, makeQuery(unique(d[[i]]$Protein), database))
        }
    }
    return(resultDf)
}

# Identical proteins translation method
.getNCBI2UniProtIP <- function(ncbiId){

    identicalProteins <- data.frame()
    tryCatch({
        identicalProteins <- getNCBIIdenticalProteins(ncbiId, format="dataframe")
    },error=function(e){})

    if(nrow(identicalProteins)>0){
        translationsRefSeq <- .getNCBI2UniProtBatch(identicalProteins[which(identicalProteins$Source=='RefSeq'),] ,"RefSeq_Protein" )
        translationsCds <- .getNCBI2UniProtBatch(identicalProteins[which(identicalProteins$Source=='INSDC'),],"EMBL-GenBank-DDBJ_CDS" )
        if(nrow(translationsRefSeq)>0 & nrow(translationsCds)>0){
            translations <- unique(rbind(translationsRefSeq, translationsCds)$Entry)
            return(c(translations))
        }else{
            return(character(0))
        }
    }else{
        return(character(0))
    }
}

# NCBI Protein to UniProt translation
getNCBIProtein2UniProt <- function(ncbiId, exhaustiveMapping = FALSE, byIdenticalProteins = TRUE){
    .checkNCBIProteinIdExists(ncbiId)
    .checkBoolean(byIdenticalProteins, 'byIdenticalProteins')
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')

    # First translation strategy
    translatedIDs <- .getNCBI2UniProtDT(ncbiId)

    # Second translation strategy
    if(byIdenticalProteins){
        translatedIDs <- c(translatedIDs, .getNCBI2UniProtIP(ncbiId))
    }

    if(!exhaustiveMapping&!is.na(unique(translatedIDs)[1])){return(unique(translatedIDs)[1])}
    else{return(unique(translatedIDs))}
}

# NCBI Nucleotide to UniProt translation
getNCBINucleotide2UniProt <- function(ncbiId, exhaustiveMapping = FALSE, byIdenticalProteins = TRUE){
    .checkNCBINucleotideIdExists(ncbiId)
    .checkBoolean(byIdenticalProteins, 'byIdenticalProteins')
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')

    # First translation strategy
    translatedIDs <- .getNCBI2UniProtDT(ncbiId)

    # Second translation strategy
    if(byIdenticalProteins){
        translatedIDs <- c(translatedIDs, .getNCBI2UniProtIP(ncbiId))
    }

    if(!exhaustiveMapping&!is.na(unique(translatedIDs)[1])){return(unique(translatedIDs)[1])}
    else{return(unique(translatedIDs))}
}

# NCBI Gene to UniProt translation
getNCBIGene2UniProt <- function(ncbiId, exhaustiveMapping = FALSE, byIdenticalProteins = TRUE){
    .checkNCBIGeneIdExists(ncbiId)
    .checkBoolean(byIdenticalProteins, 'byIdenticalProteins')
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')

    # First translation strategy
    translatedIDs <- .getNCBI2UniProtDT(ncbiId)

    # Second translation strategy
    if(byIdenticalProteins){
        translatedIDs <- c(translatedIDs, .getNCBI2UniProtIP(ncbiId))
    }

    if(!exhaustiveMapping&!is.na(unique(translatedIDs)[1])){return(unique(translatedIDs)[1])}
    else{return(unique(translatedIDs))}
}

##########################
# NCBI databases to KEGG #
##########################

# Direct translation method
.getNCBI2KEGGDT <- function(ncbiId){
    query <- keggConv("genes", paste('ncbi-proteinid:', strsplit(ncbiId,'.', fixed = T)[[1]], sep=''))

    if(!identical(query, character(0))&!identical(query, NULL)&!identical(query, NA)){
        return(c(query[[1]]))
    }else{
        return(character(0))
    }
}

# ------ Notes to facilitate documentation writing
# -- Se ha puesto desactivable lo de byIdenticalProteins y bySimilarGenes,
# Al desactivarlos, reducimos las probabilidades de traducir el gen, pero también tendrá un coste computacional menor (necesario en algunos
# casos donde por ejemplo, el método de proteínas similares tarde muchísimo)
# -- Se ha añadido un exhaustiveMapping booleano que activo, devuelve todas las posibles traducciones tanto directas como por clusters de similitud de UniProt;
# desactivada, devuelve la primera traducción directa que encuentre o, en caso contrario, la de mayor % de identidad.
# Through UniProt database method
.getNCBI2KEGGTUP <- function(ncbiId, ncbiDB, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){
    upIDs <- switch (ncbiDB,
                     'protein' = getNCBIProtein2UniProt(ncbiId, exhaustiveMapping = exhaustiveMapping, byIdenticalProteins = byIdenticalProteins),
                     'nucleotide' = getNCBINucleotide2UniProt(ncbiId, exhaustiveMapping = exhaustiveMapping, byIdenticalProteins = byIdenticalProteins),
                     'gene' = getNCBIGene2UniProt(ncbiId, exhaustiveMapping = exhaustiveMapping, byIdenticalProteins = byIdenticalProteins),
    )

    translations <- list()
    if(length(upIDs)>0){

        for(upId in upIDs){
            keggId <- getUniProt2KEGG(upId, exhaustiveMapping = exhaustiveMapping, detailedMapping = TRUE, bySimilarGenes = bySimilarGenes)
            if(length(keggId)>0){
                translations <- .mergeNamedLists(translations, keggId)
            }
        }
    }

    if(!detailedMapping){
        translations <- unique(unlist(translations, recursive = FALSE, use.names = FALSE))
        if(is.null(translations)){translations <- character(0)}
    }else{
        translations <- lapply(translations, unique)
    }
    return(translations)
}

.getNCBI2KEGG <- function(ncbiId, ncbiDB, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){
    switch(ncbiDB,
         'protein' = .checkNCBIProteinIdExists(ncbiId),
         'nucleotide' = .checkNCBINucleotideIdExists(ncbiId),
         'gene' = .checkNCBIGeneIdExists(ncbiId),
    )
    .checkBoolean(exhaustiveMapping, 'exhaustiveMapping')
    .checkBoolean(detailedMapping, 'detailedMapping')
    .checkBoolean(byIdenticalProteins, 'byIdenticalProteins')
    .checkBoolean(bySimilarGenes, 'bySimilarGenes')

    result <- list()
    translationsDT <- .getNCBI2KEGGDT(ncbiId)

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

    # Then, try to translate through UniProt database
    result <- .mergeNamedLists(result, .getNCBI2KEGGTUP(ncbiId, ncbiDB, exhaustiveMapping = exhaustiveMapping, detailedMapping = TRUE,
                                                        byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes))

    # To return a char vector of IDs if detailedMapping is false, simply extract
    # all IDs from result in a char vector and return it. Also, perform a 'unique' operation
    if(!detailedMapping){
        result <- unique(unlist(result, recursive = FALSE, use.names = FALSE))
        if(is.null(result)){result <- character(0)}
    }else{
        result <- lapply(result, unique)
    }
    return(result)
}

getNCBIProtein2KEGG <- function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){
    return(.getNCBI2KEGG(ncbiId, 'protein', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                         byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes))
}

getNCBINucleotide2KEGG <- function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){
    return(.getNCBI2KEGG(ncbiId, 'nucleotide', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                         byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes))
}

getNCBIGene2KEGG <- function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){
    return(.getNCBI2KEGG(ncbiId, 'gene', exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping,
                         byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes))
}



















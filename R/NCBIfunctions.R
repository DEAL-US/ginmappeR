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
        for(queryId in query[['links']][['gene_protein']]){
            proteinXml <- entrez_fetch(db = "protein", id = queryId, rettype = "xml")
            result <- append(result, xpathSApply(xmlParse(proteinXml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
        }
    }
    # Handle no possible translations case
    # if(length(result)==0){
    #     stop("The given NCBI Gene ID has no translation to NCBI Protein database")
    # }
    return(result)
}

getNCBIProtein2NCBIGene <- function(id){
    .checkNCBIProteinIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='protein', id=id, dbTo='gene')
    # Handle no possible translations case
    # if(length(query[['links']][['protein_gene']])==0){
    #     stop("The given NCBI Protein ID has no translation to NCBI Gene database")
    # }
    return(as.list(query[['links']][['protein_gene']]))
}

getNCBIProtein2NCBINucleotide <- function(id){
    .checkNCBIProteinIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='protein', id=id, dbTo='nucleotide')
    result <- list()
    # Handle multiple protein IDs case
    if(!identical(query[['links']][['protein_nuccore']], NULL)){
        for(queryId in query[['links']][['protein_nuccore']]){
            proteinXml <- entrez_fetch(db = "nucleotide", id = queryId, rettype = "xml")
            result <- append(result, xpathSApply(xmlParse(proteinXml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
        }
    }
    # Handle no possible translations case
    # if(length(result)==0){
    #     stop("The given NCBI Protein ID has no translation to NCBI Nucleotide database")
    # }
    return(result)
}

getNCBINucleotide2NCBIProtein <- function(id){
    .checkNCBINucleotideIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='nucleotide', id=id, dbTo='protein')
    result <- list()
    # Handle multiple protein IDs case
    if(!identical(query[['links']][['nuccore_protein']], NULL)){
        for(queryId in query[['links']][['nuccore_protein']]){
            proteinXml <- entrez_fetch(db = "protein", id = queryId, rettype = "xml")
            result <- append(result, xpathSApply(xmlParse(proteinXml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
        }
    }
    # Handle no possible translations case
    # if(length(result)==0){
    #     stop("The given NCBI Nucleotide ID has no translation to NCBI Protein database")
    # }
    return(result)
}

getNCBIGene2NCBINucleotide <- function(id){
    .checkNCBIGeneIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='gene', id=id, dbTo='nucleotide')
    result <- list()
    # Handle multiple protein IDs case
    if(!identical(query[['links']][['gene_nuccore']], NULL)){
        for(queryId in query[['links']][['gene_nuccore']]){
            proteinXml <- entrez_fetch(db = "nucleotide", id = queryId, rettype = "xml")
            result <- append(result, xpathSApply(xmlParse(proteinXml),'//GBSet/GBSeq/GBSeq_primary-accession',xmlValue)[[1]])
        }
    }
    # Handle no possible translations case
    # if(length(result)==0){
    #     stop("The given NCBI Gene ID has no translation to NCBI Nucleotide database")
    # }
    return(result)
}

getNCBINucleotide2NCBIGene <- function(id){
    .checkNCBINucleotideIdExists(id)
    query <- .getNCBIDatabasesLinks(dbFrom='nucleotide', id=id, dbTo='gene')
    # Handle no possible translations case
    # if(length(query[['links']][['nuccore_gene']])==0){
    #     stop("The given NCBI Nucleotide ID has no translation to NCBI Gene database")
    # }
    return(as.list(query[['links']][['nuccore_gene']]))
}

#############################
# NCBI databases to UniProt #
#############################

# Direct translation method
.getNCBI2UniProtDT <- function(ncbiId){
    query <- list()
    query <- append(query, suppressWarnings(mapUniProt(
        from='EMBL-GenBank-DDBJ_CDS',
        to='UniProtKB',
        query=c(ncbiId),
        columns=c('accession')
    )))
    query <- append(query, suppressWarnings(mapUniProt(
        from='RefSeq_Protein',
        to='UniProtKB',
        query=c(ncbiId),
        columns=c('accession')
    )))
    query <- append(query, suppressWarnings(mapUniProt(
        from='RefSeq_Nucleotide',
        to='UniProtKB',
        query=c(ncbiId),
        columns=c('accession')
    )))
    return(as.list(query$Entry))
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
               'ids'= return(list()),
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
               'ids'= return(as.list(identicalProteins$Protein)),
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
            return(as.list(translations))
        }else{
            return(list())
        }
    }else{
        return(list())
    }
}

# NCBI Protein to UniProt translation
getNCBIProtein2UniProt <- function(ncbiId, byIdenticalProteins = TRUE){
    .checkNCBIProteinIdExists(ncbiId)
    .checkBoolean(byIdenticalProteins, 'byIdenticalProteins')

    # First translation strategy
    translatedIDs <- .getNCBI2UniProtDT(ncbiId)

    # Second translation strategy
    if(byIdenticalProteins){
        translatedIDs <- append(translatedIDs, .getNCBI2UniProtIP(ncbiId))
    }

    # Handle no possible translations case
    # if(length(translatedIDs)==0){
    #     stop("The given NCBI Protein ID has no translation to UniProt database")
    # }else{
        return(unique(translatedIDs))
    # }
}

# NCBI Nucleotide to UniProt translation
getNCBINucleotide2UniProt <- function(ncbiId, byIdenticalProteins = TRUE){
    .checkNCBINucleotideIdExists(ncbiId)
    .checkBoolean(byIdenticalProteins, 'byIdenticalProteins')

    # First translation strategy
    translatedIDs <- .getNCBI2UniProtDT(ncbiId)

    # Second translation strategy
    if(byIdenticalProteins){
        translatedIDs <- append(translatedIDs, .getNCBI2UniProtIP(ncbiId))
    }

    # Handle no possible translations case
    # if(length(translatedIDs)==0){
        # stop("The given NCBI Nucleotide ID has no translation to UniProt database")
    # }else{
        return(unique(translatedIDs))
    # }
}

# NCBI Gene to UniProt translation
getNCBIGene2UniProt <- function(ncbiId, byIdenticalProteins = TRUE){
    .checkNCBIGeneIdExists(ncbiId)
    .checkBoolean(byIdenticalProteins, 'byIdenticalProteins')

    # First translation strategy
    translatedIDs <- .getNCBI2UniProtDT(ncbiId)

    # Second translation strategy
    if(byIdenticalProteins){
        translatedIDs <- append(translatedIDs, .getNCBI2UniProtIP(ncbiId))
    }

    # Handle no possible translations case
    # if(length(translatedIDs)==0){
    #     # stop("The given NCBI Gene ID has no translation to UniProt database")
    # }else{
    return(unique(translatedIDs))
    # }
}

##########################
# NCBI databases to KEGG #
##########################

# Direct translation method
.getNCBI2KEGGDT <- function(ncbiId){
    query <- keggConv("genes", paste('ncbi-proteinid:', strsplit(ncbiId,'.', fixed = T)[[1]], sep=''))

    if(!identical(query, character(0))&!identical(query, NULL)&!identical(query, NA)){
        return(as.list(query[[1]]))
    }else{
        return(list())
    }
}

# ------ Notes to facilitate documentation writing
# -- Se ha puesto desactivable lo de byIdenticalProteins y bySimilarGenes,
# Al desactivarlos, reducimos las probabilidades de traducir el gen, pero también tendrá un coste computacional menor (necesario en algunos
# casos donde por ejemplo, el método de proteínas similares tarde muchísimo)
# -- Se ha añadido un allTranslations booleano que activo, devuelve todas las posibles traducciones tanto directas como por clusters de similitud de UniProt;
# desactivada, devuelve la primera traducción directa que encuentre o, en caso contrario, la de mayor % de identidad.
# Through UniProt database method
.getNCBI2KEGGTUP <- function(ncbiId, ncbiDatabase, allTranslations = TRUE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){
    .getSGTFirstElement <- function(x) {
        as.numeric(x$SGT[[1]])
    }

    upIDs <- switch (ncbiDatabase,
        'protein' = getNCBIProtein2UniProt(ncbiId, byIdenticalProteins = byIdenticalProteins),
        'nucleotide' = getNCBINucleotide2UniProt(ncbiId, byIdenticalProteins = byIdenticalProteins),
        'gene' = getNCBIGene2UniProt(ncbiId, byIdenticalProteins = byIdenticalProteins),
    )

    translations <- list()
    if(length(upIDs)>0){

        # If user wants all possible translations
        if(allTranslations){
            for(upId in upIDs){
                keggId <- getUniProt2KEGG(upId, bySimilarGenes = bySimilarGenes)
                if(!identical(keggId, list('DT'=list(), 'SGT'=list()))){
                    translations[upId] <- list(keggId)
                }
            }
            if(!identical(translations, list())){
                return(translations[order(sapply(translations, .getSGTFirstElement), decreasing = TRUE)])
            }
        }else{
            # If user wants the best/faster translation
            for(upId in upIDs){
                keggId <- getUniProt2KEGG(upId, bySimilarGenes = bySimilarGenes)
                if(!identical(keggId, list('DT'=list(), 'SGT'=list()))){
                    if(length(keggId$DT)>0){
                        return(setNames(list(keggId), `upId`))
                    }
                    translations[upId] <- list(keggId)
                }
            }
            if(!identical(translations, list())){
                return(translations[order(sapply(translations, .getSGTFirstElement), decreasing = TRUE)][[1]])
            }
        }
    }
    return(translations)
}

# TODO: del método .getNCBI2KEGGTUP, en el método superior, quedarte con la DT si hay y si no, ya coger la
# de Similar Genes y devolver la explicación de 0.9, blablabla. Poner un allCoincidences o algo así, que sea
# True o False y de esta forma devolver la primera traducción o todas las posibles, según la prisa que tenga
# el usuario. Poner tmb el byIdenticalProteins y bySimilarGenes con sus correspondientes validaciones













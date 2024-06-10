library(plumber)
library(ginmappeR)

#* @apiTitle Gintegrator API
#* @apiDescription Gintegrator API provides functionalities for translating gene or protein identifiers across several biological sequence databases, including CARD, NCBI Protein, Nucleotide and Gene, UniProt, and KEGG. It allows for easy mapping of identifiers across these databases, ensuring access to the most accurate and up-to-date information. In addition to identifier mapping, Gintegrator API also offers the retrieval of NCBI identical proteins or UniProt similar genes clusters.

#* @filter cors
cors <- function(req, res) {

    res$setHeader("Access-Control-Allow-Origin", "*")

    if (req$REQUEST_METHOD == "OPTIONS") {
        res$setHeader("Access-Control-Allow-Methods","*")
        res$setHeader("Access-Control-Allow-Headers", req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS)
        res$status <- 200
        return(list())
    } else {
        plumber::forward()
    }

}


##############################################
#        Error handler helper function       #
##############################################

errorHandler <- function(func, res, ...){
    tryCatch({
        func(...)
    },error=function(error){
        if(grepl('not registered', error$message, fixed = TRUE)){
            res$status <- 404
            return(list(error='404 - Not Found', message=error$message))
        }else if(grepl('not a valid UniProt', error$message, fixed = TRUE)){
            res$status <- 400
            return(list(error='400 - Bad Request', message=error$message))
        }else{
            res$status <- 500
            return(list(error="500 - Internal Server Error", message=error$message))
        }
    })
}


###################################
#         CARD endpoints          #
###################################

#* Translate CARD to NCBI Protein
#* @tag CARD
#* @get /card/<cardId>/ncbiProtein
function(cardId, res){
    errorHandler(func = ginmappeR:::getCARD2NCBIProtein, res = res, cardId)
}

#* Translate CARD to NCBI Nucleotide
#* @tag CARD
#* @get /card/<cardId>/ncbiNucleotide
function(cardId, res){
    errorHandler(func = ginmappeR:::getCARD2NCBINucleotide, res = res, cardId)
}

#* Translate CARD to NCBI Gene
#* @tag CARD
#* @param exhaustiveMapping:bool
#* @get /card/<cardId>/ncbiGene
function(cardId, exhaustiveMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    errorHandler(func = ginmappeR:::getCARD2NCBIGene, res = res, cardId, exhaustiveMapping)
}

#* Translate CARD to UniProt
#* @tag CARD
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @get /card/<cardId>/uniprot
function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    errorHandler(func = ginmappeR:::getCARD2UniProt, res = res, cardId, exhaustiveMapping, detailedMapping)
}

#* Translate CARD to KEGG
#* @tag CARD
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param byIdenticalProteins:bool
#* @param bySimilarGenes:bool
#* @get /card/<cardId>/kegg
function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getCARD2KEGG, res = res, cardId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes)
}

###################################
#         NCBI endpoints          #
###################################

#* Translate NCBI Gene to NCBI Protein
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @get /ncbiGene/<ncbiId>/ncbiProtein
function(ncbiId, exhaustiveMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    errorHandler(func = ginmappeR:::getNCBIGene2NCBIProtein, res = res, ncbiId, exhaustiveMapping)
}

#* Translate NCBI Protein to NCBI Gene
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @get /ncbiProtein/<ncbiId>/ncbiGene
function(ncbiId, exhaustiveMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    errorHandler(func = ginmappeR:::getNCBIProtein2NCBIGene, res = res, ncbiId, exhaustiveMapping)
}

#* Translate NCBI Protein to NCBI Nucleotide
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @get /ncbiProtein/<ncbiId>/ncbiNucleotide
function(ncbiId, exhaustiveMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    errorHandler(func = ginmappeR:::getNCBIProtein2NCBINucleotide, res = res, ncbiId, exhaustiveMapping)
}

#* Translate NCBI Nucleotide to NCBI Protein
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @get /ncbiNucleotide/<ncbiId>/ncbiProtein
function(ncbiId, exhaustiveMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    errorHandler(func = ginmappeR:::getNCBINucleotide2NCBIProtein, res = res, ncbiId, exhaustiveMapping)
}

#* Translate NCBI Gene to NCBI Nucleotide
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @get /ncbiGene/<ncbiId>/ncbiNucleotide
function(ncbiId, exhaustiveMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    errorHandler(func = ginmappeR:::getNCBIGene2NCBINucleotide, res = res, ncbiId, exhaustiveMapping)
}

#* Translate NCBI Nucleotide to NCBI Gene
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @get /ncbiNucleotide/<ncbiId>/ncbiGene
function(ncbiId, exhaustiveMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    errorHandler(func = ginmappeR:::getNCBINucleotide2NCBIGene, res = res, ncbiId, exhaustiveMapping)
}

#* Get NCBI identical proteins
#* @tag NCBI
#* @param format
#* @get /identicalProteins/<ncbiId>
function(ncbiId, format = 'ids'){
    if (!format %in% c("ids", "dataframe")){format <- 'ids'}
    ginmappeR:::getNCBIIdenticalProteins(ncbiId, format)
}

#* Translate NCBI Protein to UniProt
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param byIdenticalProteins:bool
#* @get /ncbiProtein/<ncbiId>/uniprot
function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    errorHandler(func = ginmappeR:::getNCBIProtein2UniProt, res = res, ncbiId, exhaustiveMapping, detailedMapping, byIdenticalProteins)
}

#* Translate NCBI Nucleotide to UniProt
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param byIdenticalProteins:bool
#* @get /ncbiNucleotide/<ncbiId>/uniprot
function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    errorHandler(func = ginmappeR:::getNCBINucleotide2UniProt, res = res, ncbiId, exhaustiveMapping, detailedMapping, byIdenticalProteins)
}

#* Translate NCBI Gene to UniProt
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param byIdenticalProteins:bool
#* @get /ncbiGene/<ncbiId>/uniprot
function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    errorHandler(func = ginmappeR:::getNCBIGene2UniProt, res = res, ncbiId, exhaustiveMapping, detailedMapping, byIdenticalProteins)
}

#* Translate NCBI Protein to KEGG
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param byIdenticalProteins:bool
#* @param bySimilarGenes:bool
#* @get /ncbiProtein/<ncbiId>/kegg
function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getNCBIProtein2KEGG, res = res, ncbiId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes)
}

#* Translate NCBI Nucleotide to KEGG
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param byIdenticalProteins:bool
#* @param bySimilarGenes:bool
#* @get /ncbiNucleotide/<ncbiId>/kegg
function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getNCBINucleotide2KEGG, res = res, ncbiId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes)
}

#* Translate NCBI Gene to KEGG
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param byIdenticalProteins:bool
#* @param bySimilarGenes:bool
#* @get /ncbiGene/<ncbiId>/kegg
function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getNCBIGene2KEGG, res = res, ncbiId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes)
}


#* Translate NCBI Protein to CARD
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @get /ncbiProtein/<ncbiId>/card
function(ncbiId, exhaustiveMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    errorHandler(func = ginmappeR:::getNCBIProtein2CARD, res = res, ncbiId, exhaustiveMapping)
}

#* Translate NCBI Nucleotide to CARD
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @get /ncbiNucleotide/<ncbiId>/card
function(ncbiId, exhaustiveMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    errorHandler(func = ginmappeR:::getNCBINucleotide2CARD, res = res, ncbiId, exhaustiveMapping)
}

#* Translate NCBI Gene to CARD
#* @tag NCBI
#* @param exhaustiveMapping:bool
#* @get /ncbiGene/<ncbiId>/card
function(ncbiId, exhaustiveMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    errorHandler(func = ginmappeR:::getNCBIGene2CARD, res = res, ncbiId, exhaustiveMapping)
}

###################################
#        UniProt endpoints        #
###################################

#* Get UniProt similar genes
#* @tag UniProt
#* @param clusterIdentity
#* @param clusterNames:bool
#* @get /similarGenes/<upId>
function(upId, clusterIdentity = '1.0', clusterNames = FALSE){
    if (!clusterIdentity %in% c("1.0", "0.9", "0.5")){clusterIdentity <- '1.0'}
    if (clusterNames %in% c("TRUE", "FALSE")){clusterNames <- as.logical(clusterNames)} else {clusterNames <- FALSE}
    ginmappeR:::getUniProtSimilarGenes(upId, clusterIdentity, clusterNames)
}

#* Translate UniProt to KEGG
#* @tag UniProt
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /uniprot/<upId>/kegg
function(upId, exhaustiveMapping = FALSE, bySimilarGenes = TRUE, detailedMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getUniProt2KEGG, res = res, upId, exhaustiveMapping, bySimilarGenes, detailedMapping)
}

#* Translate UniProt to NCBI Protein
#* @tag UniProt
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /uniprot/<upId>/ncbiProtein
function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getUniProt2NCBIProtein, res = res, upId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate UniProt to NCBI Nucleotide
#* @tag UniProt
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /uniprot/<upId>/ncbiNucleotide
function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getUniProt2NCBINucleotide, res = res, upId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate UniProt to NCBI Gene
#* @tag UniProt
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /uniprot/<upId>/ncbiGene
function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getUniProt2NCBIGene, res = res, upId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate UniProt to CARD
#* @tag UniProt
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /uniprot/<upId>/card
function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getUniProt2CARD, res = res, upId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}


###################################
#         KEGG endpoints          #
###################################

#* Translate KEGG to UniProt
#* @tag KEGG
#* @param exhaustiveMapping:bool
#* @get /kegg/<keggId>/uniprot
function(keggId, exhaustiveMapping = FALSE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    errorHandler(func = ginmappeR:::getKEGG2UniProt, res = res, keggId, exhaustiveMapping)
}

#* Translate KEGG to NCBI Protein
#* @tag KEGG
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /kegg/<keggId>/ncbiProtein
function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getKEGG2NCBIProtein, res = res, keggId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate KEGG to NCBI Nucleotide
#* @tag KEGG
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /kegg/<keggId>/ncbiNucleotide
function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getKEGG2NCBINucleotide, res = res, keggId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate KEGG to NCBI Gene
#* @tag KEGG
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /kegg/<keggId>/ncbiGene
function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getKEGG2NCBIGene, res = res, keggId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate KEGG to CARD
#* @tag KEGG
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /kegg/<keggId>/card
function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE, res){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    errorHandler(func = ginmappeR:::getKEGG2CARD, res = res, keggId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}




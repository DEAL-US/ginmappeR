library(plumber)
library(ginmappeR)

#* @apiTitle ginmappeR API

###################################
#        Utility endpoints        #
###################################

#* Updates CARD database version
#* @get /updateCard
ginmappeR:::updateCARDDataBase

###################################
#         CARD endpoints          #
###################################

#* Translate CARD to NCBI Protein
#* @get /card/<cardId>/ncbiProtein
ginmappeR:::getCARD2NCBIProtein

#* Translate CARD to NCBI Nucleotide
#* @get /card/<cardId>/ncbiNucleotide
ginmappeR:::getCARD2NCBINucleotide

#* Translate CARD to NCBI Gene
#* @param exhaustiveMapping:bool
#* @get /card/<cardId>/ncbiGene
function(cardId, exhaustiveMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    ginmappeR:::getCARD2NCBIGene(cardId, exhaustiveMapping)
}

#* Translate CARD to UniProt
#* @param exhaustiveMapping:bool
#* @get /card/<cardId>/uniprot
function(cardId, exhaustiveMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    ginmappeR:::getCARD2UniProt(cardId, exhaustiveMapping)
}

#* Translate CARD to KEGG
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param byIdenticalProteins:bool
#* @param bySimilarGenes:bool
#* @get /card/<cardId>/kegg
function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getCARD2KEGG(cardId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes)
}

###################################
#         KEGG endpoints          #
###################################

#* Translate KEGG to UniProt
#* @get /kegg/<keggId>/uniprot
ginmappeR:::getKEGG2UniProt

#* Translate KEGG to NCBI Protein
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /kegg/<keggId>/ncbiProtein
function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getKEGG2NCBIProtein(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate KEGG to NCBI Nucleotide
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /kegg/<keggId>/ncbiNucleotide
function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getKEGG2NCBINucleotide(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate KEGG to NCBI Gene
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /kegg/<keggId>/ncbiGene
function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getKEGG2NCBIGene(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate KEGG to CARD
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /kegg/<keggId>/card
function(keggId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getKEGG2CARD(keggId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

###################################
#         NCBI endpoints          #
###################################

#* Translate NCBI Gene to NCBI Protein
#* @param exhaustiveMapping:bool
#* @get /ncbiGene/<ncbiId>/ncbiProtein
function(ncbiId, exhaustiveMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    ginmappeR:::getNCBIGene2NCBIProtein(ncbiId, exhaustiveMapping)
}

#* Translate NCBI Protein to NCBI Gene
#* @param exhaustiveMapping:bool
#* @get /ncbiProtein/<ncbiId>/ncbiGene
function(ncbiId, exhaustiveMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    ginmappeR:::getNCBIProtein2NCBIGene(ncbiId, exhaustiveMapping)
}

#* Translate NCBI Protein to NCBI Nucleotide
#* @param exhaustiveMapping:bool
#* @get /ncbiProtein/<ncbiId>/ncbiNucleotide
function(ncbiId, exhaustiveMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    ginmappeR:::getNCBIProtein2NCBINucleotide(ncbiId, exhaustiveMapping)
}

#* Translate NCBI Nucleotide to NCBI Protein
#* @param exhaustiveMapping:bool
#* @get /ncbiNucleotide/<ncbiId>/ncbiProtein
function(ncbiId, exhaustiveMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    ginmappeR:::getNCBINucleotide2NCBIProtein(ncbiId, exhaustiveMapping)
}

#* Translate NCBI Gene to NCBI Nucleotide
#* @param exhaustiveMapping:bool
#* @get /ncbiGene/<ncbiId>/ncbiNucleotide
function(ncbiId, exhaustiveMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    ginmappeR:::getNCBIGene2NCBINucleotide(ncbiId, exhaustiveMapping)
}

#* Translate NCBI Nucleotide to NCBI Gene
#* @param exhaustiveMapping:bool
#* @get /ncbiNucleotide/<ncbiId>/ncbiGene
function(ncbiId, exhaustiveMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    ginmappeR:::getNCBINucleotide2NCBIGene(ncbiId, exhaustiveMapping)
}

#* Get NCBI identical proteins
#* @param format
#* @get /identicalProteins/<ncbiId>
function(ncbiId, format = 'ids'){
    if (!format %in% c("ids", "dataframe")){format <- 'ids'}
    ginmappeR:::getNCBIIdenticalProteins(ncbiId, format)
}

#* Translate NCBI Protein to UniProt
#* @param exhaustiveMapping:bool
#* @param byIdenticalProteins:bool
#* @get /ncbiProtein/<ncbiId>/uniprot
function(ncbiId, exhaustiveMapping = FALSE, byIdenticalProteins = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    ginmappeR:::getNCBIProtein2UniProt(ncbiId, exhaustiveMapping, byIdenticalProteins)
}

#* Translate NCBI Nucleotide to UniProt
#* @param exhaustiveMapping:bool
#* @param byIdenticalProteins:bool
#* @get /ncbiNucleotide/<ncbiId>/uniprot
function(ncbiId, exhaustiveMapping = FALSE, byIdenticalProteins = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    ginmappeR:::getNCBINucleotide2UniProt(ncbiId, exhaustiveMapping, byIdenticalProteins)
}

#* Translate NCBI Gene to UniProt
#* @param exhaustiveMapping:bool
#* @param byIdenticalProteins:bool
#* @get /ncbiGene/<ncbiId>/uniprot
function(ncbiId, exhaustiveMapping = FALSE, byIdenticalProteins = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    ginmappeR:::getNCBIGene2UniProt(ncbiId, exhaustiveMapping, byIdenticalProteins)
}

#* Translate NCBI Protein to KEGG
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param byIdenticalProteins:bool
#* @param bySimilarGenes:bool
#* @get /ncbiProtein/<ncbiId>/kegg
function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getNCBIProtein2KEGG(ncbiId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes)
}

#* Translate NCBI Nucleotide to KEGG
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param byIdenticalProteins:bool
#* @param bySimilarGenes:bool
#* @get /ncbiNucleotide/<ncbiId>/kegg
function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getNCBINucleotide2KEGG(ncbiId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes)
}

#* Translate NCBI Gene to KEGG
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param byIdenticalProteins:bool
#* @param bySimilarGenes:bool
#* @get /ncbiGene/<ncbiId>/kegg
function(ncbiId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (byIdenticalProteins %in% c("TRUE", "FALSE")){byIdenticalProteins <- as.logical(byIdenticalProteins)} else {byIdenticalProteins <- TRUE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getNCBIGene2KEGG(ncbiId, exhaustiveMapping, detailedMapping, byIdenticalProteins, bySimilarGenes)
}


#* Translate NCBI Protein to CARD
#* @param exhaustiveMapping:bool
#* @get /ncbiProtein/<ncbiId>/card
function(ncbiId, exhaustiveMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    ginmappeR:::getNCBIProtein2CARD(ncbiId, exhaustiveMapping)
}

#* Translate NCBI Nucleotide to CARD
#* @param exhaustiveMapping:bool
#* @get /ncbiNucleotide/<ncbiId>/card
function(ncbiId, exhaustiveMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    ginmappeR:::getNCBINucleotide2CARD(ncbiId, exhaustiveMapping)
}

#* Translate NCBI Gene to CARD
#* @param exhaustiveMapping:bool
#* @get /ncbiGene/<ncbiId>/card
function(ncbiId, exhaustiveMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    ginmappeR:::getNCBIGene2CARD(ncbiId, exhaustiveMapping)
}

###################################
#        UniProt endpoints        #
###################################

#* Get UniProt similar genes
#* @param clusterIdentity
#* @param clusterNames:bool
#* @get /similarGenes/<upId>
function(upId, clusterIdentity = '1.0', clusterNames = FALSE){
    if (!clusterIdentity %in% c("1.0", "0.9", "0.5")){clusterIdentity <- '1.0'}
    if (clusterNames %in% c("TRUE", "FALSE")){clusterNames <- as.logical(clusterNames)} else {clusterNames <- FALSE}
    ginmappeR:::getUniProtSimilarGenes(upId, clusterIdentity, clusterNames)
}

#* Translate UniProt to KEGG
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /uniprot/<upId>/kegg
function(upId, exhaustiveMapping = FALSE, bySimilarGenes = TRUE, detailedMapping = FALSE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getUniProt2KEGG(upId, exhaustiveMapping, bySimilarGenes, detailedMapping)
}

#* Translate UniProt to NCBI Protein
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /uniprot/<upId>/ncbiProtein
function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getUniProt2NCBIProtein(upId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate UniProt to NCBI Nucleotide
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /uniprot/<upId>/ncbiNucleotide
function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getUniProt2NCBINucleotide(upId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate UniProt to NCBI Gene
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /uniprot/<upId>/ncbiGene
function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getUniProt2NCBIGene(upId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}

#* Translate UniProt to CARD
#* @param exhaustiveMapping:bool
#* @param detailedMapping:bool
#* @param bySimilarGenes:bool
#* @get /uniprot/<upId>/card
function(upId, exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = TRUE){
    if (exhaustiveMapping %in% c("TRUE", "FALSE")){exhaustiveMapping <- as.logical(exhaustiveMapping)} else {exhaustiveMapping <- FALSE}
    if (detailedMapping %in% c("TRUE", "FALSE")){detailedMapping <- as.logical(detailedMapping)} else {detailedMapping <- FALSE}
    if (bySimilarGenes %in% c("TRUE", "FALSE")){bySimilarGenes <- as.logical(bySimilarGenes)} else {bySimilarGenes <- TRUE}
    ginmappeR:::getUniProt2CARD(upId, exhaustiveMapping, detailedMapping, bySimilarGenes)
}





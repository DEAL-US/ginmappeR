
###################################
# CARD database to NCBI databases #
###################################

getCARD2NCBIProtein <- function(cardId){

    .checkIfCARDIsDownloaded()
    .checkCARDIdExists(cardId)

    aroIndex <- read.csv(paste(get_cache_dir(get_pkg_info("ginmappeR")),'/card-data/aro_index.tsv', sep=''), sep='\t')
    proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession
    return(proteinId)
}

getCARD2NCBINucleotide <- function(cardId){

    .checkIfCARDIsDownloaded()
    .checkCARDIdExists(cardId)

    aroIndex <- read.csv(paste(get_cache_dir(get_pkg_info("ginmappeR")),'/card-data/aro_index.tsv', sep=''), sep='\t')
    nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
    return(nucleotideId)
}

getCARD2NCBIGene <- function(cardId, exhaustiveMapping = FALSE){

    .checkIfCARDIsDownloaded()
    .checkCARDIdExists(cardId)
    .checkBoolean(exhaustiveMapping)

    aroIndex <- read.csv(paste(get_cache_dir(get_pkg_info("ginmappeR")),'/card-data/aro_index.tsv', sep=''), sep='\t')
    nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
    proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession

    result <- unique(c(getNCBIProtein2NCBIGene(proteinId, exhaustiveMapping = exhaustiveMapping), getNCBINucleotide2NCBIGene(nucleotideId, exhaustiveMapping = exhaustiveMapping)))

    if(!exhaustiveMapping & length(result)>0){
        return(result[1])
    }
    return(result)
}

############################
# CARD database to UniProt #
############################

getCARD2UniProt <- function(cardId, exhaustiveMapping = FALSE){

    .checkIfCARDIsDownloaded()
    .checkCARDIdExists(cardId)

    aroIndex <- read.csv(paste(get_cache_dir(get_pkg_info("ginmappeR")),'/card-data/aro_index.tsv', sep=''), sep='\t')
    nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
    proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession

    result <- getNCBIProtein2UniProt(proteinId, exhaustiveMapping = exhaustiveMapping, byIdenticalProteins = TRUE)
    if(exhaustiveMapping | length(result)==0){
        result <- unique(c(result, getNCBINucleotide2UniProt(nucleotideId, exhaustiveMapping = exhaustiveMapping, byIdenticalProteins = TRUE)))
    }
    return(result)
}

#########################
# CARD database to KEGG #
#########################

### TODO: documentation:
# Warn that byIdenticalProteins and bySimilarGenes combined with exhaustiveMapping takes a long time to finish
getCARD2KEGG <- function(cardId, exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = TRUE, bySimilarGenes = TRUE){

    .checkIfCARDIsDownloaded()
    .checkCARDIdExists(cardId)

    aroIndex <- read.csv(paste(get_cache_dir(get_pkg_info("ginmappeR")),'/card-data/aro_index.tsv', sep=''), sep='\t')
    nucleotideId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$DNA.Accession
    proteinId <- aroIndex[aroIndex$ARO.Accession == (if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)),]$Protein.Accession

    result <- getNCBIProtein2KEGG(proteinId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes)

    if(exhaustiveMapping | length(result)==0){
        aux <- getNCBINucleotide2KEGG(nucleotideId, exhaustiveMapping = exhaustiveMapping, detailedMapping = detailedMapping, byIdenticalProteins = byIdenticalProteins, bySimilarGenes = bySimilarGenes)
        if(detailedMapping){
            result <- lapply(.mergeNamedLists(result, aux), unique)
        }else{
            result <- unique(c(result, aux))
        }
    }
    return(result)
}





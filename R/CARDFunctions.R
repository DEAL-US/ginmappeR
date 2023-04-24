
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







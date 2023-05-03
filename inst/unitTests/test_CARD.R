library(RUnit)
library(ginmappeR)
library('UniProt.ws')
library('KEGGREST')
library('httr')
library('rentrez')
library('XML')
library('pkgfilecache')

source('../R/CARDFunctions.R')
source('../R/utilsFunctions.R')

print('-------------')
print(getwd())
print('-------------')

#########################
# CARD database to NCBI #
#########################

### Test getCARD2NCBIProtein

# Positive cases
checkEquals(getCARD2NCBIProtein('3002535'), 'CAA38525.1')
checkEquals(getCARD2NCBIProtein('ARO:3002535'), 'CAA38525.1')
checkEquals(getCARD2NCBIProtein('3003988'), 'APB03221.1')
# Incorrect CARD ID case
checkException(getCARD2NCBIProtein('test'))

### Test getCARD2NCBINucleotide

# Positive cases
checkEquals(getCARD2NCBINucleotide('3002535'), 'X54723.1')
checkEquals(getCARD2NCBINucleotide('ARO:3002535'), 'X54723.1')
checkEquals(getCARD2NCBINucleotide('3003988'), 'KX531051.1')
# Incorrect CARD ID case
checkException(getCARD2NCBINucleotide('test'))

### Test getCARD2NCBIGene

# Positive cases
checkEquals(getCARD2NCBIGene('3002524'), c('29426913'))
checkEquals(getCARD2NCBIGene('ARO:3002524'), c('29426913'))
checkEquals(getCARD2NCBIGene('ARO:3002524', TRUE), c('29426913'))
checkEquals(getCARD2NCBIGene('3002525'), c('886648'))
# No translation case
checkEquals(getCARD2NCBIGene('3005061'), character(0))

############################
# CARD database to UniProt #
############################

### Test getCARD2UniProt

# Positive cases
checkEquals(getCARD2UniProt('3002867'), c('Q9ZIF9'))
checkTrue(length(getCARD2UniProt('3002867', TRUE)) == 38)
checkEquals(getCARD2UniProt('3003649'), c('A0A0K0TQH5'))
# No translation case
checkEquals(getCARD2UniProt('3006267'), character(0))

#########################
# CARD database to KEGG #
#########################

### Test getCARD2KEGG

# Positive cases
checkEquals(getCARD2KEGG('3000938'), c('ag:AAF19151'))
checkEquals(getCARD2KEGG('3000938', detailedMapping = TRUE), list('DT'=c('ag:AAF19151')))
checkEquals(getCARD2KEGG('3001109', detailedMapping = TRUE), list('0.9'=c('ag:BAA84973')))
checkEquals(getCARD2KEGG('3002511', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE),
            list('0.5'=c("chk:D4L85_28045","proe:H9L23_08075")))
# No translation cases
checkEquals(getCARD2KEGG('3006267', detailedMapping = FALSE, bySimilarGenes = FALSE), character(0))
checkEquals(getCARD2KEGG('3006267', detailedMapping = TRUE, bySimilarGenes = FALSE), list())








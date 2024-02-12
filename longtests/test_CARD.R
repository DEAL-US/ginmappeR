library(RUnit)
library(ginmappeR)
library('UniProt.ws')
library('KEGGREST')
library('httr')
library('rentrez')
library('XML')

utils::globalVariables('cardPath')
cardPath <<- tempdir()

# GitHub Actions imports
# source('../../../R/CARDFunctions.R')
# source('../../../R/utilsFunctions.R')

# # Local execution imports
setwd('../00_pkg_src/ginmappeR/')
source('R/CARDFunctions.R')
source('R/utilsFunctions.R')

#########################
# CARD database to NCBI #
#########################

### Test getCARD2NCBIProtein
message('Testing getCARD2NCBIProtein')
# Positive cases
.testEquals(getCARD2NCBIProtein('3002535'), 'CAA38525.1')
.testEquals(getCARD2NCBIProtein('ARO:3002535'), 'CAA38525.1')
.testEquals(getCARD2NCBIProtein('3003988'), 'APB03221.1')
.testEquals(getCARD2NCBIProtein(c('3003988', 'ARO:3002535')), c('APB03221.1', 'CAA38525.1'))
# Incorrect CARD ID case
.testEquals(getCARD2NCBIProtein('test'), NA)
.testEquals(getCARD2NCBIProtein(c('test', 'test')), c(NA, NA))


### Test getCARD2NCBINucleotide
message('Testing getCARD2NCBINucleotide')
# Positive cases
.testEquals(getCARD2NCBINucleotide('3002535'), 'X54723.1')
.testEquals(getCARD2NCBINucleotide('ARO:3002535'), 'X54723.1')
.testEquals(getCARD2NCBINucleotide('3003988'), 'KX531051.1')
.testEquals(getCARD2NCBINucleotide(c('3003988', 'ARO:3002535')), c('KX531051.1', 'X54723.1'))
# Incorrect CARD ID case
.testEquals(getCARD2NCBINucleotide('test'), NA)
.testEquals(getCARD2NCBINucleotide(c('test','test')), c(NA,NA))


### Test getCARD2NCBIGene
message('Testing getCARD2NCBIGene')
# Positive cases
.testEquals(getCARD2NCBIGene('3002525'), c('886648'))
.testEquals(getCARD2NCBIGene(c('3002525','3002525')), c('886648', '886648'))
# # No translation case
# .testEquals(getCARD2NCBIGene('3005061'), NA)
# Incorrect CARD ID case
.testEquals(getCARD2NCBIGene('test'), NA)

############################
# CARD database to UniProt #
############################

### Test getCARD2UniProt
message('Testing getCARD2UniProt')
# Positive cases
.testEquals(getCARD2UniProt('3002867'), c('Q9ZIF9'))
.testEquals(getCARD2UniProt('3002867', detailedMapping=TRUE), list('DT'=c('Q9ZIF9')))
.testEquals(getCARD2UniProt(c('3002867', '3002867'), detailedMapping=TRUE), c('DT'=c('Q9ZIF9'),'DT'=c('Q9ZIF9')))
# checkTrue(length(getCARD2UniProt('3002867', TRUE)) == 35)
# .testEquals(getCARD2UniProt('3003649'), c('A0A0K0TQH5'))
# No translation case
.testEquals(getCARD2UniProt('3006267'), NA)
.testEquals(getCARD2UniProt('3006267', detailedMapping=TRUE), NA)
# Incorrect CARD ID case
.testEquals(getCARD2UniProt('test'), NA)

#########################
# CARD database to KEGG #
#########################

### Test getCARD2KEGG
message('Testing getCARD2KEGG')
# Positive cases
.testEquals(getCARD2KEGG('3000938'), c('ag:AAF19151'))
.testEquals(getCARD2KEGG('3000938', detailedMapping = TRUE), list('DT'='ag:AAF19151'))
.testEquals(getCARD2KEGG(c('3000938','3000938'), detailedMapping = TRUE), c('DT'=c('ag:AAF19151'), 'DT'=c('ag:AAF19151')))
# .testEquals(getCARD2KEGG('3001109', detailedMapping = TRUE), list('0.9'=c('ag:BAA84973')))
# .testEquals(getCARD2KEGG('3002511', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE),
            # list('0.5'=c("chk:D4L85_28045","proe:H9L23_08075"))) # Takes a long time
# No translation cases
.testEquals(getCARD2KEGG('3006267', detailedMapping = FALSE, bySimilarGenes = FALSE), NA)
.testEquals(getCARD2KEGG('3006267', detailedMapping = TRUE, bySimilarGenes = FALSE), NA)
# Incorrect CARD ID case
.testEquals(getCARD2KEGG('test'), NA)







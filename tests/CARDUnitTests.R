library(RUnit)
source('../R/CARDFunctions.R')
source('../R/utilsFunctions.R')
library('UniProt.ws')
library('KEGGREST')
library('httr')
library('rentrez')
library('XML')

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
checkEquals(getCARD2NCBIGene('3002524'), '29426913')
checkEquals(getCARD2NCBIGene('ARO:3002524'), '29426913')
checkEquals(getCARD2NCBIGene('ARO:3002524', TRUE), c('29426913'))
checkEquals(getCARD2NCBIGene('3002525'), '886648')
# No translation case
checkEquals(getCARD2NCBIGene('3005061'), character(0))


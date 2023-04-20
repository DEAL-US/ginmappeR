library(RUnit)
source('../R/KEGGFunctions.R')
source('../R/utilsFunctions.R')
library('UniProt.ws')
library('KEGGREST')
library('httr')
library('rentrez')
library('XML')

###################################
# KEGG database auxiliar function #
###################################

### Test .checkKEGGIdExists

# Positive case
checkEquals(.checkKEGGIdExists('lacy:A4V08_06425'), NULL)
# ID not registered case
checkException(.checkKEGGIdExists('example_id'))
# ID not valid case
checkException(.checkKEGGIdExists('lacy:A4V08_'))


############################
# KEGG database to UniProt #
############################

### Test getKEGG2UniProt

# Positive cases
checkEquals(getKEGG2UniProt('lacy:A4V08_06425'), c('A0A1C7G297'))
checkEquals(getKEGG2UniProt('mct:MCR_1292'), c('D5VCZ9'))
checkEquals(getKEGG2UniProt('abc:ACICU_00223'), c('A0A7U3XVN2'))
# No translation case
checkEquals(getKEGG2UniProt('eco:b0204'), character(0))

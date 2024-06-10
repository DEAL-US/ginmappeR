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
# source('../../../R/KEGGFunctions.R')
# source('../../../R/00utils.R')

# # Local execution imports

# source('../R/KEGGFunctions.R')
# source('../R/00utils.R')

###################################
# KEGG database auxiliar function #
###################################

### Test .checkKEGGIdExists
# message('Testing .checkKEGGIdExists')
# # Positive case
# .testEquals(.checkKEGGIdExists('lacy:A4V08_06425'), NULL)
# # ID not registered case
# RUnit::checkException(.checkKEGGIdExists('example_id'))
# # ID not valid case
# RUnit::checkException(.checkKEGGIdExists('lacy:A4V08_'))

############################
# KEGG database to UniProt #
############################

### Test getKEGG2UniProt
message('Testing getKEGG2UniProt')
# Positive cases
# .testEquals(getKEGG2UniProt('abc:ACICU_00223'), c('A0A7U3XVN2'))
ginmappeR:::.testEquals(getKEGG2UniProt('aag:5579347'), c('A0A1S4G4Z1'))
ginmappeR:::.testEquals(getKEGG2UniProt('aag:5579347', exhaustiveMapping = TRUE), list(c('A0A1S4G4Z1', 'Q1DGP2')))
ginmappeR:::.testEquals(getKEGG2UniProt(c('aag:5579347','aag:5579347')), c('A0A1S4G4Z1','A0A1S4G4Z1'))
# No translation case
ginmappeR:::.testEquals(getKEGG2UniProt('eco:b0204'), NA)
# ID not valid case
# ginmappeR:::.testEquals(getKEGG2UniProt('test'), NA)

###################################
# KEGG database to NCBI databases #
###################################

### Test .getKEGG2NCBIDT
# message('Testing .getKEGG2NCBIDT')
# # Positive cases
# .testEquals(.getKEGG2NCBIDT('lacy:A4V08_06425'), c('ANU45512'))
# # .testEquals(.getKEGG2NCBIDT('mct:MCR_1292'), c('ADG61552'))
# # .testEquals(.getKEGG2NCBIDT('abc:ACICU_00223'), c('ACC55535'))
# .testEquals(.getKEGG2NCBIDT('aag:5579347'), c('5579347','XP_001647653'))

### Test .getKEGG2NCBITUP
# message('Testing .getKEGG2NCBITUP')
# Positive cases
# Sys.sleep(3)
# .testEquals(.getKEGG2NCBITUP('aag:5579347', 'protein', bySimilarGenes = FALSE), list('DT'=c("EAT32321.1","XP_001647653.1")))
# .testEquals(.getKEGG2NCBITUP('aag:5579347', 'nucleotide', bySimilarGenes = FALSE), list())
# Sys.sleep(3)
# .testEquals(.getKEGG2NCBITUP('aag:5579347', 'gene', bySimilarGenes = FALSE), list('DT'=c("EAT32321.1","XP_001647653.1")))
# .testEquals(.getKEGG2NCBITUP('aag:5579347', 'nucleotide', bySimilarGenes = TRUE), list('0.5'='XM_001230804.1'))
# .testEquals(.getKEGG2NCBITUP('llo:LLO_2673', 'nucleotide', exhaustiveMapping = TRUE, bySimilarGenes = TRUE),
#             list('DT'=c('NC_013861.1'), '0.5'=c("NZ_CCVW01000004.1","NZ_LNYO01000024.1","NZ_UGNZ01000001.1","NZ_UASS01000015.1",
#                                                 "NZ_UGOX01000001.1","NZ_UGGV01000001.1"))) # Takes a long time
# .testEquals(.getKEGG2NCBITUP('llo:LLO_2673', 'gene', bySimilarGenes = TRUE), list('DT'=c('CBJ13102.1')))

### Test getKEGG2NCBIProtein, getKEGG2NCBINucleotide, getKEGG2NCBIGene
message('Testing getKEGG2NCBIProtein, getKEGG2NCBINucleotide, getKEGG2NCBIGene')
# Positive cases
# checkTrue(length(getKEGG2NCBIProtein('llo:LLO_2673', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE)$`0.5`)==51)
# .testEquals(getKEGG2NCBINucleotide('llo:LLO_2673', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE),
#             list('DT'=c('NC_013861.1'), '0.5'=c("NZ_CCVW01000004.1","NZ_LNYO01000024.1","NZ_UGNZ01000001.1","NZ_UASS01000015.1",
#                                                 "NZ_UGOX01000001.1","NZ_UGGV01000001.1"))) # Takes a long time
ginmappeR:::.testEquals(getKEGG2NCBINucleotide('llo:LLO_2673', bySimilarGenes = TRUE), c('NC_013861.1'))
ginmappeR:::.testEquals(getKEGG2NCBINucleotide(c('llo:LLO_2673','llo:LLO_2673'), bySimilarGenes = TRUE), c('NC_013861.1', 'NC_013861.1'))
# .testEquals(getKEGG2NCBINucleotide('aag:5579347', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE), list('0.5'='XM_001230804.1'))
# .testEquals(getKEGG2NCBIGene('aag:5579347', detailedMapping = TRUE, bySimilarGenes = TRUE), list('DT'='5579347'))
ginmappeR:::.testEquals(getKEGG2NCBIGene('aag:5579347', bySimilarGenes = TRUE), c('5579347'))
ginmappeR:::.testEquals(getKEGG2NCBIGene(c('aag:5579347','aag:5579347'), bySimilarGenes = TRUE), c('5579347', '5579347'))
# ID not valid case
# ginmappeR:::.testEquals(getKEGG2NCBIProtein('test'), NA)
# ginmappeR:::.testEquals(getKEGG2NCBINucleotide('test'), NA)
# ginmappeR:::.testEquals(getKEGG2NCBIGene('test'), NA)

#########################
# KEGG database to CARD #
#########################

### Test getKEGG2CARD
message('Testing getKEGG2CARD')
# Positive cases
# testEquals(getKEGG2CARD('llo:LLO_2673', detailedMapping = TRUE), list('0.5' = c('ARO:3004591'))) # Takes too long
ginmappeR:::.testEquals(getKEGG2CARD('ag:ACC85616'), c("ARO:3002804"))
ginmappeR:::.testEquals(getKEGG2CARD(c('ag:ACC85616','ag:ACC85616')), c("ARO:3002804", "ARO:3002804"))
ginmappeR:::.testEquals(getKEGG2CARD('ag:ACC85616', detailedMapping = TRUE), list('DT' = c('ARO:3002804')))
# .testEquals(getKEGG2CARD('ag:CAJ47134'), c("ARO:3001133"))
# .testEquals(getKEGG2CARD('ag:CAJ47134', detailedMapping = TRUE), list('DT' = c('ARO:3001133')))
# .testEquals(getKEGG2CARD('ag:ACC85616', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = FALSE), list('DT' = c('ARO:3002804')))
# .testEquals(getKEGG2CARD('ag:ACC85616', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE), list('DT' = c('ARO:3002804')))
# .testEquals(getKEGG2CARD('llo:LLO_2673', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE), list('0.5' = c("ARO:3004591", "ARO:3004584", "ARO:3004587", "ARO:3004586", "ARO:3004613", "ARO:3004590", "ARO:3004589", "ARO:3004582", "ARO:3004581")))
# No translation case
ginmappeR:::.testEquals(getKEGG2CARD('llo:LLO_2673', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = FALSE), list(NULL))
ginmappeR:::.testEquals(getKEGG2CARD('llo:LLO_2673', exhaustiveMapping = TRUE, detailedMapping = FALSE, bySimilarGenes = FALSE), list(NULL))
ginmappeR:::.testEquals(getKEGG2CARD('llo:LLO_2673', exhaustiveMapping = FALSE, detailedMapping = FALSE, bySimilarGenes = FALSE), NA)
ginmappeR:::.testEquals(getKEGG2CARD('llo:LLO_2673', exhaustiveMapping = FALSE, detailedMapping = TRUE, bySimilarGenes = FALSE), NA)
# ID not valid case
# ginmappeR:::.testEquals(getKEGG2CARD('test'), NA)


library(RUnit)
library(ginmappeR)
library('UniProt.ws')
library('KEGGREST')
library('httr')
library('rentrez')
library('XML')
library('pkgfilecache')
source('../../R/KEGGFunctions.R', chdir = TRUE)
source('../../R/utilsFunctions.R', chdir = TRUE)

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
checkEquals(getKEGG2UniProt('aag:5579347'), c('A0A1S4G4Z1','Q1DGP2'))
# No translation case
checkEquals(getKEGG2UniProt('eco:b0204'), character(0))

###################################
# KEGG database to NCBI databases #
###################################

### Test .getKEGG2NCBIDT

# Positive cases
checkEquals(.getKEGG2NCBIDT('lacy:A4V08_06425'), c('ANU45512'))
checkEquals(.getKEGG2NCBIDT('mct:MCR_1292'), c('ADG61552'))
checkEquals(.getKEGG2NCBIDT('abc:ACICU_00223'), c('ACC55535'))
checkEquals(.getKEGG2NCBIDT('aag:5579347'), c('5579347','XP_001647653'))

### Test .getKEGG2NCBITUP

# Positive cases
checkEquals(.getKEGG2NCBITUP('aag:5579347', 'protein', bySimilarGenes = FALSE), list('DT'=c("EAT32321.1","XP_001647653.1")))
checkEquals(.getKEGG2NCBITUP('aag:5579347', 'nucleotide', bySimilarGenes = FALSE), list())
checkEquals(.getKEGG2NCBITUP('aag:5579347', 'gene', bySimilarGenes = FALSE), list('DT'=c("EAT32321.1","XP_001647653.1")))
checkEquals(.getKEGG2NCBITUP('aag:5579347', 'nucleotide', bySimilarGenes = TRUE), list('0.5'='XM_001230804.1'))
checkEquals(.getKEGG2NCBITUP('llo:LLO_2673', 'nucleotide', exhaustiveMapping = TRUE, bySimilarGenes = TRUE),
            list('DT'=c('NC_013861.1'), '0.5'=c("NZ_CCVW01000004.1","NZ_LNYO01000024.1","NZ_UGNZ01000001.1","NZ_UASS01000015.1",
                                                "NZ_UGOX01000001.1","NZ_UGGV01000001.1")))
checkEquals(.getKEGG2NCBITUP('llo:LLO_2673', 'gene', bySimilarGenes = TRUE), list('DT'=c('CBJ13102.1')))

### Test getKEGG2NCBIProtein, getKEGG2NCBINucleotide, getKEGG2NCBIGene

# Positive cases
checkTrue(length(getKEGG2NCBIProtein('llo:LLO_2673', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE)$`0.5`)==51)
checkEquals(getKEGG2NCBINucleotide('llo:LLO_2673', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE),
            list('DT'=c('NC_013861.1'), '0.5'=c("NZ_CCVW01000004.1","NZ_LNYO01000024.1","NZ_UGNZ01000001.1","NZ_UASS01000015.1",
                                                "NZ_UGOX01000001.1","NZ_UGGV01000001.1")))
checkEquals(getKEGG2NCBINucleotide('llo:LLO_2673', bySimilarGenes = TRUE), c('NC_013861.1'))
checkEquals(getKEGG2NCBINucleotide('aag:5579347', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE), list('0.5'='XM_001230804.1'))
checkEquals(getKEGG2NCBIGene('aag:5579347', detailedMapping = TRUE, bySimilarGenes = TRUE), list('DT'='5579347'))
checkEquals(getKEGG2NCBIGene('aag:5579347', bySimilarGenes = TRUE), c('5579347'))

#########################
# KEGG database to CARD #
#########################

### Test getKEGG2CARD

# Positive cases
checkEquals(getKEGG2CARD('llo:LLO_2673', detailedMapping = TRUE), list('0.5' = c('ARO:3004591')))
checkEquals(getKEGG2CARD('ag:ACC85616'), c("ARO:3002804"))
checkEquals(getKEGG2CARD('ag:ACC85616', detailedMapping = TRUE), list('DT' = c('ARO:3002804')))
checkEquals(getKEGG2CARD('ag:CAJ47134'), c("ARO:3001133"))
checkEquals(getKEGG2CARD('ag:CAJ47134', detailedMapping = TRUE), list('DT' = c('ARO:3001133')))
checkEquals(getKEGG2CARD('ag:ACC85616', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = FALSE), list('DT' = c('ARO:3002804')))
checkEquals(getKEGG2CARD('ag:ACC85616', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE), list('DT' = c('ARO:3002804')))
checkEquals(getKEGG2CARD('llo:LLO_2673', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = TRUE), list('0.5' = c("ARO:3004591", "ARO:3004584", "ARO:3004587", "ARO:3004586", "ARO:3004613", "ARO:3004590", "ARO:3004589", "ARO:3004582", "ARO:3004581")))
# No translation case
checkEquals(getKEGG2CARD('llo:LLO_2673', exhaustiveMapping = TRUE, detailedMapping = TRUE, bySimilarGenes = FALSE), list())



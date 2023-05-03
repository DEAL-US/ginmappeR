library(RUnit)
library(ginmappeR)
library('UniProt.ws')
library('KEGGREST')
library('httr')
library('rentrez')
library('XML')
library('pkgfilecache')
# source('../R/UniProtFunctions.R')
# source('../R/utilsFunctions.R')

######################################
# UniProt database auxiliar function #
######################################

### Test .checkUniProtIdExists

# Positive case
checkEquals(.checkUniProtIdExists('P0DTH6'), NULL)
# ID not registered case
checkException(.checkUniProtIdExists('P0ZTH2'))
# ID not valid case
checkException(.checkUniProtIdExists('P0ZTUH2'))

############################
# UniProt database to KEGG #
############################

### Test .getUniProt2KEGGDT

# Positive cases
checkEquals(.getUniProt2KEGGDT('G9JVE6'), c('ag:AEX08599'))
checkEquals(.getUniProt2KEGGDT('P37639'), c('eco:b3516','ecj:JW3484'))
# No translation case
checkEquals(.getUniProt2KEGGDT('P0DTH6'), character(0))

### Test getUniProtSimilarGenes

# Positive cases
checkTrue(length(getUniProtSimilarGenes('G0L217', clusterIdentity = '0.5', clusterNames = TRUE)$UniRef50_A0A5S3PNY9)==4)
checkTrue(length(getUniProtSimilarGenes('G9JVE6', clusterIdentity = '1.0', clusterNames = TRUE)$UniRef100_G9JVE6)==10)
checkTrue(length(getUniProtSimilarGenes('G9JVE6', clusterIdentity = '1.0'))==10)
# No similar genes case
checkEquals(getUniProtSimilarGenes('G0L217', clusterIdentity = '1.0', clusterNames = TRUE), list('UniRef100_G0L217'=character(0)))
checkEquals(getUniProtSimilarGenes('G0L217', clusterIdentity = '1.0'), character(0))
# Invalid cluster identity provided
checkException(getUniProtSimilarGenes('G9JVE6','test'))

### Test .getUniProt2KEGGSGT
# Positive cases
checkEquals(.getUniProt2KEGGSGT('B2ZPD3'), list('1.0' = c('kpb:FH42_26825')))
checkEquals(.getUniProt2KEGGSGT('A0A0B5ECY2'), list('0.9' = c('ag:CAZ39946')))
checkEquals(.getUniProt2KEGGSGT('A0A2R4PHC7'), list('0.5' = c('ffl:HYN86_12595')))
checkEquals(.getUniProt2KEGGSGT('A0A2R4PHC7', TRUE), list('0.5' = c('ffl:HYN86_12595', 'fcr:HYN56_14615', 'fls:GLV81_10715')))
# No translation case
checkEquals(.getUniProt2KEGGSGT('A0A1S7BGS4'), list())

### Test getUniProt2KEGG

checkEquals(getUniProt2KEGG('G9JVE6', detailedMapping = TRUE), list('DT'=c('ag:AEX08599')))
checkEquals(getUniProt2KEGG('G9JVE6'), c('ag:AEX08599'))
checkEquals(getUniProt2KEGG('A0A2R4PHC7', detailedMapping = TRUE), list('0.5'=c('ffl:HYN86_12595')))
checkEquals(getUniProt2KEGG('A0A2R4PHC7'), c('ffl:HYN86_12595'))
checkEquals(getUniProt2KEGG('A0A2R4PHC7', exhaustiveMapping = TRUE, detailedMapping = TRUE), list('0.5'=c('ffl:HYN86_12595', 'fcr:HYN56_14615', 'fls:GLV81_10715')))
checkTrue(length(getUniProt2KEGG('Q6XL56', exhaustiveMapping = TRUE, detailedMapping = TRUE)[['0.5']])==20)
checkEquals(getUniProt2KEGG('G9JVE6', exhaustiveMapping = TRUE, bySimilarGenes = FALSE, detailedMapping = TRUE), list('DT'=c('ag:AEX08599')))
checkEquals(getUniProt2KEGG('Q6XL56', exhaustiveMapping = TRUE, bySimilarGenes = FALSE, detailedMapping = TRUE), list('DT'=c('ag:AAP69916')))
checkEquals(getUniProt2KEGG('Q6XL56', exhaustiveMapping = TRUE, bySimilarGenes = FALSE), c('ag:AAP69916'))
# No translation cases
checkEquals(getUniProt2KEGG('A0A1S7BGS4'), character(0))
checkEquals(getUniProt2KEGG('A0A1S7BGS4', detailedMapping = TRUE), list())
checkEquals(getUniProt2KEGG('A0A1S7BGS4', FALSE, FALSE), character(0))
# ID not valid case
checkException(getUniProt2KEGG('P0ZTUH2'))

############################
# UniProt database to NCBI #
############################

### Test .getNCBI2UniProtDT

# Positive cases
checkEquals(.getUniProt2NCBIDT('A0SNL9'), c("ABH10964.1","WP_002449123.1","NZ_VDQI01000026.1"))
checkEquals(.getUniProt2NCBIDT('A0SNL9 B0BKZ8'), c("ABH10964.1","CAP69659.1","WP_002449123.1","NZ_VDQI01000026.1"))
# ID not registered case
checkEquals(.getUniProt2NCBIDT('test'), character(0))

### Test .getUniProt2NCBISGT

# Positive cases
test <- .getUniProt2NCBISGT('B0BKZ8', ncbiDB = 'nucleotide')
checkEquals(.getUniProt2NCBISGT('A0SNL9', ncbiDB = 'protein'), list('0.9' = c("ACX34100.1","WP_032492560.1")))
checkEquals(.getUniProt2NCBISGT('B0BKZ8', ncbiDB = 'nucleotide'), list('1.0' = c("NZ_UWXB01000002.1")))
checkTrue(length(.getUniProt2NCBISGT('A0A6H2TXZ6', ncbiDB = 'protein', exhaustiveMapping = TRUE)$`0.5`)==15)
# No translation case
checkEquals(.getUniProt2NCBISGT('Q6R7P5', ncbiDB = 'protein'), list())

# Test getUniProt2NCBIProtein, getUniProt2NCBINucleotide, getUniProt2NCBIGene

# Positive cases
checkTrue(length(getUniProt2NCBIProtein('A0A6H2TXZ6', exhaustiveMapping = TRUE, detailedMapping = TRUE)$`0.5`)==15)
checkTrue(length(getUniProt2NCBINucleotide('A0A6H2TXZ6', exhaustiveMapping = TRUE, detailedMapping = TRUE)$`0.5`)==2)
checkEquals(getUniProt2NCBINucleotide('A0A6H2TXZ6', exhaustiveMapping = FALSE, detailedMapping = TRUE), list('0.5' = c("NZ_LMFS01000003.1")))
checkEquals(getUniProt2NCBIProtein('A0SNL9', exhaustiveMapping = FALSE, detailedMapping = FALSE), c('ABH10964.1'))
checkEquals(getUniProt2NCBINucleotide('A0SNL9', exhaustiveMapping = FALSE, detailedMapping = FALSE), c('NZ_VDQI01000026.1'))
checkEquals(getUniProt2NCBIGene('A0SNL9', exhaustiveMapping = FALSE, detailedMapping = FALSE), c('AFH57403.1'))
# No translation cases
checkEquals(getUniProt2NCBIGene('A0A1S7BGS4', bySimilarGenes = FALSE), character(0))
checkEquals(getUniProt2NCBIGene('A0A1S7BGS4', bySimilarGenes = FALSE, detailedMapping = TRUE), list())
# ID not valid case
checkException(getUniProt2NCBIProtein('P0ZTUH2'))

############################
# UniProt database to CARD #
############################

### Test getUniProt2CARD

# Positive cases
checkEquals(getUniProt2CARD('A0A1S7BGS4'), c("ARO:3004185"))
checkEquals(getUniProt2CARD('A0A4R5XW64', bySimilarGenes = TRUE), c("ARO:3004185"))
checkEquals(getUniProt2CARD('Q8GNY5', detailedMapping = TRUE), list('DT'=c("ARO:3003552")))
checkEquals(getUniProt2CARD('Q8GNY5'), c('ARO:3003552'))
checkEquals(getUniProt2CARD('A0A6H2TXZ6', exhaustiveMapping = TRUE, detailedMapping = TRUE), list('DT' = c('ARO:3005012'), '0.9'= c('ARO:3005013'), '0.5'=c('ARO:3005013')))
checkEquals(getUniProt2CARD('A0A6H2TXZ6', exhaustiveMapping = FALSE, detailedMapping = TRUE), list('DT' = c('ARO:3005012')))
# No translation cases
checkEquals(getUniProt2CARD('A0A4R5XW64', bySimilarGenes = FALSE), character(0))
checkEquals(getUniProt2CARD('A0A4R5XW64', detailedMapping = TRUE, bySimilarGenes = FALSE), list())


library(RUnit)
source('../R/UniProtFunctions.R')
source('../R/utilsFunctions.R')
library('UniProt.ws')
library('KEGGREST')
library('httr')


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




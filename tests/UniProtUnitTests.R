library(RUnit)
source('../R/UniProtFunctions.R')
source('../R/utilsFunctions.R')
library('UniProt.ws')

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
checkEquals(.getUniProt2KEGGDT('G9JVE6'), list('ag:AEX08599'))
checkEquals(.getUniProt2KEGGDT('P37639'), list('eco:b3516','ecj:JW3484'))
# No translation case
checkEquals(.getUniProt2KEGGDT('P0DTH6'), list())

### Test getUniProtSimilarGenes

# Positive cases
checkTrue(length(getUniProtSimilarGenes('G0L217', clusterIdentity = '0.5')$UniRef50_A0A5S3PNY9)==4)
checkTrue(length(getUniProtSimilarGenes('G9JVE6', clusterIdentity = '1.0')$UniRef100_G9JVE6)==10)
# No similar genes case
checkEquals(getUniProtSimilarGenes('G0L217', clusterIdentity = '1.0'), list('UniRef100_G0L217'=list()))
# Invalid cluster identity provided
checkException(getUniProtSimilarGenes('G9JVE6','test'))

### Test .getUniProt2KEGGSGT
# Positive cases
checkEquals(.getUniProt2KEGGSGT('B2ZPD3'), list('1.0','B9V4W1','kpb:FH42_26825'))
checkEquals(.getUniProt2KEGGSGT('A0A0B5ECY2'), list('0.9','C7C422','ag:CAZ39946'))
checkEquals(.getUniProt2KEGGSGT('A0A2R4PHC7'), list('0.5','A0A344LTZ2','ffl:HYN86_12595'))
# No translation case
checkEquals(.getUniProt2KEGGSGT('A0A1S7BGS4'), list())

### Test getUniProt2KEGG

# Positive cases
checkEquals(getUniProt2KEGG('G9JVE6'), list( 'DT' = list('ag:AEX08599'), 'SGT' = list('0.9', 'C7C422', 'ag:CAZ39946')))
checkEquals(getUniProt2KEGG('G9JVE6', FALSE), list('DT'=list('ag:AEX08599')))
# No translation cases
checkEquals(getUniProt2KEGG('A0A1S7BGS4'), list('DT'=list(), 'SGT'=list()))
checkEquals(getUniProt2KEGG('A0A1S7BGS4', FALSE), list('DT'=list()))
# ID not valid case
checkException(getUniProt2KEGG('P0ZTUH2'))






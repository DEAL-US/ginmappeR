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

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
# source('../../../R/UniProtFunctions.R')
# source('../../../R/utilsFunctions.R')

# # Local execution imports
# setwd('../00_pkg_src/ginmappeR/')
# source('../00_pkg_src/ginmappeR/R/UniProtFunctions.R')
# source('../00_pkg_src/ginmappeR/R/utilsFunctions.R')

######################################
# UniProt database auxiliar function #
######################################

# ### Test .checkUniProtIdExists
# message('Testing .checkUniProtIdExists')
# # Positive case
# ..testEquals(.checkUniProtIdExists('P0DTH6'), NULL)
# # ID not registered case
# checkException(.checkUniProtIdExists('P0ZTH2'))
# # ID not valid case
# checkException(.checkUniProtIdExists('P0ZTUH2'))

############################
# UniProt database to KEGG #
############################

### Test .getUniProt2KEGGDT
# message('Testing .getUniProt2KEGGDT')
# Positive cases
# .testEquals(.getUniProt2KEGGDT('G9JVE6'), c('ag:AEX08599'))
# .testEquals(.getUniProt2KEGGDT('P37639'), c('eco:b3516','ecj:JW3484'))
# No translation case
# .testEquals(.getUniProt2KEGGDT('P0DTH6'), character(0))

### Test getUniProtSimilarGenes
message('Testing getUniProtSimilarGenes')
# Positive cases
RUnit::checkTrue(length(getUniProtSimilarGenes('G0L217', clusterIdentity = '0.5', clusterNames = TRUE)[[1]]$UniRef50_A0A6L9EH79)==3)
# checkTrue(length(getUniProtSimilarGenes('G9JVE6', clusterIdentity = '1.0', clusterNames = TRUE)$UniRef100_G9JVE6)==10)
# checkTrue(length(getUniProtSimilarGenes('G9JVE6', clusterIdentity = '1.0'))==10)
# No similar genes case
ginmappeR:::.testEquals(getUniProtSimilarGenes('G0L217', clusterIdentity = '1.0', clusterNames = TRUE), list(list('UniRef100_G0L217'=NULL)))
ginmappeR:::.testEquals(getUniProtSimilarGenes(c('G0L217','G0L217'), clusterIdentity = '1.0', clusterNames = TRUE), list(list('UniRef100_G0L217'=NULL),list('UniRef100_G0L217'=NULL)))
ginmappeR:::.testEquals(getUniProtSimilarGenes(c('G0L217','G0L217'), clusterIdentity = '0.5', clusterNames = TRUE),
                        list(list('UniRef50_A0A6L9EH79'=c("A0A6L9EH79","A0A967B2L0","A0A7X3D206")), list('UniRef50_A0A6L9EH79'=c("A0A6L9EH79","A0A967B2L0","A0A7X3D206"))))
ginmappeR:::.testEquals(getUniProtSimilarGenes('G0L217', clusterIdentity = '1.0'), list(NULL))
# Invalid cluster identity provided
RUnit::checkException(getUniProtSimilarGenes('G9JVE6','test'))
# ID not valid case
ginmappeR:::.testEquals(getUniProtSimilarGenes('test'), list(NULL))
ginmappeR:::.testEquals(getUniProtSimilarGenes(c('test', 'test')), list(NULL, NULL))


### Test .getUniProt2KEGGSGT
# message('Testing .getUniProt2KEGGSGT')
# Positive cases
# .testEquals(.getUniProt2KEGGSGT('B2ZPD3'), list('1.0' = c('kpb:FH42_26825')))
# .testEquals(.getUniProt2KEGGSGT('A0A0B5ECY2'), list('0.9' = c('ag:CAZ39946')))
# .testEquals(.getUniProt2KEGGSGT('A0A2R4PHC7'), list('0.5' = c('fcr:HYN56_14615')))
# .testEquals(.getUniProt2KEGGSGT('A0A2R4PHC7', TRUE), list('0.5' = c('ffl:HYN86_12595', 'fcr:HYN56_14615', 'fls:GLV81_10715')))

### Test getUniProt2KEGG
message('Testing getUniProt2KEGG')
# ginmappeR:::.testEquals(getUniProt2KEGG('G9JVE6', detailedMapping = TRUE), list('DT'=c('ag:AEX08599')))
# ginmappeR:::.testEquals(getUniProt2KEGG('G9JVE6'), c('ag:AEX08599'))
# ginmappeR:::.testEquals(getUniProt2KEGG(c('G9JVE6','G9JVE6')), list(c('ag:AEX08599'),c('ag:AEX08599')))
# .testEquals(getUniProt2KEGG('A0A2R4PHC7', detailedMapping = TRUE), list('0.5'=c('fcr:HYN56_14615')))
# .testEquals(getUniProt2KEGG('A0A2R4PHC7'), c('fcr:HYN56_14615'))
# .testEquals(getUniProt2KEGG('A0A2R4PHC7', exhaustiveMapping = TRUE, detailedMapping = TRUE), list('0.5'=c('ffl:HYN86_12595', 'fcr:HYN56_14615', 'fls:GLV81_10715')))
# checkTrue(length(getUniProt2KEGG('Q6XL56', exhaustiveMapping = TRUE, detailedMapping = TRUE)[['0.5']])==20)
# .testEquals(getUniProt2KEGG('G9JVE6', exhaustiveMapping = TRUE, bySimilarGenes = FALSE, detailedMapping = TRUE), list('DT'=c('ag:AEX08599')))
# .testEquals(getUniProt2KEGG('Q6XL56', exhaustiveMapping = TRUE, bySimilarGenes = FALSE, detailedMapping = TRUE), list('DT'=c('ag:AAP69916')))
# .testEquals(getUniProt2KEGG('Q6XL56', exhaustiveMapping = TRUE, bySimilarGenes = FALSE), c('ag:AAP69916'))

# ginmappeR:::.testEquals(getUniProt2KEGG('P0ZTUH2'), NA)

# ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5')), NA)
# ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5'), exhaustiveMapping = TRUE), list(NA))
# ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5'), detailedMapping = TRUE, exhaustiveMapping = FALSE), NA)
# ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5'), detailedMapping = TRUE, exhaustiveMapping = TRUE), list(NA))
#
# ginmappeR:::.testEquals(getUniProt2KEGG(c('test')), NA)
# ginmappeR:::.testEquals(getUniProt2KEGG(c('test'), detailedMapping = TRUE), NA)
# ginmappeR:::.testEquals(getUniProt2KEGG(c('test'), exhaustiveMapping = TRUE), list(NA))
# ginmappeR:::.testEquals(getUniProt2KEGG(c('test'), detailedMapping = TRUE, exhaustiveMapping = TRUE), list(NA))
#
# ginmappeR:::.testEquals(getUniProt2KEGG(c('A0A2R4PHC7')), c("fcs:TRV642_2906"))
# ginmappeR:::.testEquals(getUniProt2KEGG(c('A0A2R4PHC7'), exhaustiveMapping = TRUE), list(c("fcs:TRV642_2906", "fcr:HYN56_14615",
#                                                                                     "fls:GLV81_10715", "ffl:HYN86_12595")))
# ginmappeR:::.testEquals(getUniProt2KEGG(c('A0A2R4PHC7'), detailedMapping = TRUE), list('0.5'="fcs:TRV642_2906"))
# ginmappeR:::.testEquals(getUniProt2KEGG(c('A0A2R4PHC7'), detailedMapping = TRUE, exhaustiveMapping = TRUE),
#                        list(list('0.5'=c("fcs:TRV642_2906", "fcr:HYN56_14615","fls:GLV81_10715", "ffl:HYN86_12595"))))
#
# ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5', 'test','A0A2R4PHC7')), c(NA,NA, "fcs:TRV642_2906"))
# ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5', 'test','A0A2R4PHC7'), exhaustiveMapping = TRUE), list(NA, NA,
#                                                  c("fcs:TRV642_2906", "fcr:HYN56_14615","fls:GLV81_10715", "ffl:HYN86_12595")))
# ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5', 'test','A0A2R4PHC7'), detailedMapping = TRUE), c(NA,NA,'0.5'="fcs:TRV642_2906"))
# ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5', 'test','A0A2R4PHC7'), detailedMapping = TRUE, exhaustiveMapping = TRUE),
#                        list(NA,NA,list('0.5'= c("fcs:TRV642_2906", "fcr:HYN56_14615","fls:GLV81_10715", "ffl:HYN86_12595"))))


ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5')), NA)
ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5'), exhaustiveMapping = TRUE), list(NULL))
ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5'), detailedMapping = TRUE, exhaustiveMapping = FALSE), NA)
ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5'), detailedMapping = TRUE, exhaustiveMapping = TRUE), list(NULL))

ginmappeR:::.testEquals(getUniProt2KEGG(c('test')), NA)
ginmappeR:::.testEquals(getUniProt2KEGG(c('test'), detailedMapping = TRUE), NA)
ginmappeR:::.testEquals(getUniProt2KEGG(c('test'), exhaustiveMapping = TRUE), list(NULL))
ginmappeR:::.testEquals(getUniProt2KEGG(c('test'), detailedMapping = TRUE, exhaustiveMapping = TRUE), list(NULL))

ginmappeR:::.testEquals(getUniProt2KEGG(c(character(0))), character(0))
ginmappeR:::.testEquals(getUniProt2KEGG(c(character(0)), detailedMapping = TRUE), character(0))
ginmappeR:::.testEquals(getUniProt2KEGG(c(character(0)), exhaustiveMapping = TRUE), list())
ginmappeR:::.testEquals(getUniProt2KEGG(c(character(0)), detailedMapping = TRUE, exhaustiveMapping = TRUE), list())

ginmappeR:::.testEquals(getUniProt2KEGG(c('A0A2R4PHC7')), c("fcs:TRV642_2906"))
ginmappeR:::.testEquals(getUniProt2KEGG(c('A0A2R4PHC7'), exhaustiveMapping = TRUE), list(c("fcs:TRV642_2906", "fcr:HYN56_14615",
                                                                                           "fls:GLV81_10715", "ffl:HYN86_12595")))
ginmappeR:::.testEquals(getUniProt2KEGG(c('A0A2R4PHC7'), detailedMapping = TRUE), list('0.5'="fcs:TRV642_2906"))
ginmappeR:::.testEquals(getUniProt2KEGG(c('A0A2R4PHC7'), detailedMapping = TRUE, exhaustiveMapping = TRUE),
                        list(list('0.5'=c("fcs:TRV642_2906", "fcr:HYN56_14615","fls:GLV81_10715", "ffl:HYN86_12595"))))

ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5', 'test', character(0), 'A0A2R4PHC7')), c(NA, NA, "fcs:TRV642_2906"))
ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5', 'test', character(0), 'A0A2R4PHC7'), detailedMapping = TRUE), c(NA,NA,'0.5'="fcs:TRV642_2906"))
ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5', 'test', character(0),'A0A2R4PHC7'), exhaustiveMapping = TRUE), list(NULL, NULL,
                                                                                c("fcs:TRV642_2906", "fcr:HYN56_14615","fls:GLV81_10715", "ffl:HYN86_12595")))
ginmappeR:::.testEquals(getUniProt2KEGG(c('Q6R7P5', 'test', character(0), 'A0A2R4PHC7', character(0)), detailedMapping = TRUE, exhaustiveMapping = TRUE),
                        list(NULL,NULL, list('0.5'= c("fcs:TRV642_2906", "fcr:HYN56_14615","fls:GLV81_10715", "ffl:HYN86_12595"))))

############################
# UniProt database to NCBI #
############################

### Test .getNCBI2UniProtDT
# message('Testing .getNCBI2UniProtDT')
# Positive cases
# .testEquals(.getUniProt2NCBIDT('A0SNL9'), c("ABH10964.1","WP_002449123.1","NZ_VDQI01000026.1"))
# .testEquals(.getUniProt2NCBIDT('A0SNL9 B0BKZ8'), c("ABH10964.1","CAP69659.1","WP_002449123.1","NZ_VDQI01000026.1"))
# ID not registered case
# .testEquals(.getUniProt2NCBIDT('test'), character(0))

### Test .getUniProt2NCBISGT
# message('Testing .getUniProt2NCBISGT')
# Positive cases
# .testEquals(.getUniProt2NCBISGT('A0SNL9', ncbiDB = 'protein'), list('0.9' = c("ACX34100.1","WP_032492560.1")))
# .testEquals(.getUniProt2NCBISGT('B0BKZ8', ncbiDB = 'nucleotide'), list('1.0' = c("NZ_UWXB01000002.1")))
# checkTrue(length(.getUniProt2NCBISGT('A0A6H2TXZ6', ncbiDB = 'protein', exhaustiveMapping = TRUE)$`0.5`)==15)
# No translation case
# .testEquals(.getUniProt2NCBISGT('Q6R7P5', ncbiDB = 'protein'), list())

# Test getUniProt2NCBIProtein, getUniProt2NCBINucleotide, getUniProt2NCBIGene
message('Testing getUniProt2NCBIProtein, getUniProt2NCBINucleotide, getUniProt2NCBIGene')
# Positive cases

ginmappeR:::.testEquals(length(getUniProt2NCBIProtein('A0A6H2TXZ6', exhaustiveMapping = TRUE, detailedMapping = TRUE)[[1]]$`0.5`), 10)
ginmappeR:::.testEquals(length(getUniProt2NCBINucleotide('A0A6H2TXZ6', exhaustiveMapping = TRUE, detailedMapping = TRUE)[[1]]$`0.5`), 1)
ginmappeR:::.testEquals(length(getUniProt2NCBINucleotide(c('A0A6H2TXZ6', 'A0A6H2TXZ6'), exhaustiveMapping = TRUE, detailedMapping = TRUE)[[1]]$`0.5`), 1)
ginmappeR:::.testEquals(getUniProt2NCBINucleotide(c('A0A6H2TXZ6', 'A0A6H2TXZ6'), exhaustiveMapping = FALSE, detailedMapping = FALSE), c('NZ_CP013692.1','NZ_CP013692.1'))
# ginmappeR:::.testEquals(getUniProt2NCBIProtein('A0SNL9', exhaustiveMapping = FALSE, detailedMapping = FALSE), c('ABH10964.1'))
# ginmappeR:::.testEquals(getUniProt2NCBIProtein(c('A0SNL9','A0SNL9'), exhaustiveMapping = FALSE, detailedMapping = FALSE), list(c('ABH10964.1'),c('ABH10964.1')))
# ginmappeR:::.testEquals(getUniProt2NCBINucleotide('A0SNL9', exhaustiveMapping = FALSE, detailedMapping = FALSE), c('NZ_VDQI01000026.1'))
# ginmappeR:::.testEquals(getUniProt2NCBIGene('A0SNL9', exhaustiveMapping = FALSE, detailedMapping = FALSE), c('AFH57403.1'))
# No translation cases
ginmappeR:::.testEquals(getUniProt2NCBIGene('A0A1S7BGS4', bySimilarGenes = FALSE), NA)
# .testEquals(getUniProt2NCBIGene('A0A1S7BGS4', bySimilarGenes = FALSE, detailedMapping = TRUE), list())
# ID not valid case
# ginmappeR:::.testEquals(getUniProt2NCBIProtein('test'), NA)
# ginmappeR:::.testEquals(getUniProt2NCBINucleotide('test'), NA)
# ginmappeR:::.testEquals(getUniProt2NCBIGene('test'), NA)

############################
# UniProt database to CARD #
############################

### Test getUniProt2CARD
message('Testing getUniProt2CARD')
# Positive cases
ginmappeR:::.testEquals(getUniProt2CARD('A0A1S7BGS4'), c("ARO:3004185"))
ginmappeR:::.testEquals(getUniProt2CARD(c('A0A1S7BGS4','A0A1S7BGS4')), c("ARO:3004185", "ARO:3004185"))
# .testEquals(getUniProt2CARD('A0A4R5XW64', bySimilarGenes = TRUE), c("ARO:3004185"))
ginmappeR:::.testEquals(getUniProt2CARD(c('Q8GNY5', 'Q8GNY5'), detailedMapping = TRUE), c('DT'=c("ARO:3003552"), 'DT'=c("ARO:3003552")))
# .testEquals(getUniProt2CARD('Q8GNY5'), c('ARO:3003552'))
# .testEquals(getUniProt2CARD('A0A6H2TXZ6', exhaustiveMapping = TRUE, detailedMapping = TRUE), list('DT' = c('ARO:3005012'), '0.9'= c('ARO:3005013'), '0.5'=c('ARO:3005013')))
ginmappeR:::.testEquals(getUniProt2CARD('A0A6H2TXZ6', exhaustiveMapping = FALSE, detailedMapping = TRUE), list('DT' = c('ARO:3005012')))
# No translation cases
ginmappeR:::.testEquals(getUniProt2CARD('A0A4R5XW64', bySimilarGenes = FALSE), NA)
# .testEquals(getUniProt2CARD('A0A4R5XW64', detailedMapping = TRUE, bySimilarGenes = FALSE), list())
# ID not valid case
# ginmappeR:::.testEquals(getUniProt2CARD('test'), NA)

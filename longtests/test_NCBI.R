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
# source('../../../R/NCBIFunctions.R')
# source('../../../R/UniProtFunctions.R')
# source('../../../R/utilsFunctions.R')

# # Local execution imports
setwd('../00_pkg_src/ginmappeR/')
source('R/NCBIFunctions.R')
source('R/UniProtFunctions.R')
source('R/utilsFunctions.R')

#####################################
# NCBI databases inter-translations #
#####################################

### Test getNCBIGene2NCBIProtein
message('Testing getNCBIGene2NCBIProtein')
# Positive case
.testEquals(getNCBIGene2NCBIProtein('76524190'), c('WP_001082319'))
.testEquals(getNCBIGene2NCBIProtein(c('76524190', '76524190')), c('WP_001082319','WP_001082319'))
.testEquals(getNCBIGene2NCBIProtein('76524190', exhaustiveMapping = TRUE), list(c('WP_001082319')))
# ID not registered case
# .testEquals(getNCBIGene2NCBIProtein('test'), NA)
# No translation case
# .testEquals(getNCBIGene2NCBIProtein('WP_001082319'), NA)

### Test getNCBIProtein2NCBIGene
message('Testing getNCBIProtein2NCBIGene')
# Positive case
.testEquals(getNCBIProtein2NCBIGene('CAA79696'), c('1272'))
.testEquals(getNCBIProtein2NCBIGene(c('CAA79696','CAA79696')), c('1272', '1272'))
.testEquals(getNCBIProtein2NCBIGene('CAA79696', exhaustiveMapping = TRUE), list(c('1272')))
.testEquals(getNCBIProtein2NCBIGene(c('CAA79696','CAA79696'), exhaustiveMapping = TRUE), list(c('1272'), c('1272')))
# ID not registered case
# .testEquals(getNCBIProtein2NCBIGene('test'), NA)
# No translation case
# .testEquals(getNCBIProtein2NCBIGene('WP_011997479'), NA)

### Test getNCBIProtein2NCBINucleotide
message('Testing getNCBIProtein2NCBINucleotide')
# Positive case
.testEquals(getNCBIProtein2NCBINucleotide('AFH35853'), c('JQ394987'))
.testEquals(getNCBIProtein2NCBINucleotide(c('AFH35853','AFH35853')), c('JQ394987', 'JQ394987'))
.testEquals(getNCBIProtein2NCBINucleotide('AFH35853', exhaustiveMapping = TRUE), list(c('JQ394987')))
# ID not registered case
# .testEquals(getNCBIProtein2NCBINucleotide('test'), NA)

### Test getNCBINucleotide2NCBIProtein
message('Testing getNCBINucleotide2NCBIProtein')
# Positive case
.testEquals(getNCBINucleotide2NCBIProtein('JQ394987'), c('AFH35853'))
.testEquals(getNCBINucleotide2NCBIProtein(c('JQ394987','JQ394987')), c('AFH35853', 'AFH35853'))
# ID not registered case
# .testEquals(getNCBINucleotide2NCBIProtein('test'), NA)

### Test getNCBIGene2NCBINucleotide
message('Testing getNCBIGene2NCBINucleotide')
# Positive case
.testEquals(getNCBIGene2NCBINucleotide('76524190'), c('NZ_CP059690'))
.testEquals(getNCBIGene2NCBINucleotide(c('76524190','76524190')), c('NZ_CP059690', 'NZ_CP059690'))
.testEquals(getNCBIGene2NCBINucleotide('76524190', exhaustiveMapping = TRUE), list(c('NZ_CP059690')))
# ID not registered case
# .testEquals(getNCBIGene2NCBINucleotide('test'), NA)
# No translation case
# .testEquals(getNCBIGene2NCBINucleotide('WP_001082319'), NA)

### Test getNCBINucleotide2NCBIGene
message('Testing getNCBINucleotide2NCBIGene')
# Positive case
.testEquals(getNCBINucleotide2NCBIGene('Z21488'), c('1272'))
.testEquals(getNCBINucleotide2NCBIGene(c('Z21488','Z21488')), c('1272', '1272'))
.testEquals(getNCBINucleotide2NCBIGene('Z21488', exhaustiveMapping = TRUE), list(c('1272')))
# ID not registered case
# .testEquals(getNCBINucleotide2NCBIGene('test'), NA)
# No translation case
# .testEquals(getNCBINucleotide2NCBIGene('KF513177'), NA)

#############################
# NCBI databases to UniProt #
#############################

### Test .getNCBI2UniProtDT
# message('Testing .getNCBI2UniProtDT')
# # Positive cases
# .testEquals(.getNCBI2UniProtDT('AEJ08681'), c('F8TCS6'))
# .testEquals(.getNCBI2UniProtDT('AEJ08681 AGQ48857.1 CAA79696'), c('F8TCS6','S5FUH0','Q12860'))
# # .testEquals(.getNCBI2UniProtDT('CAA79696'), c('Q12860'))
# # ID not registered case
# .testEquals(.getNCBI2UniProtDT('test'), character(0))

### Test getNCBIIdenticalProteins
message('Testing getNCBIIdenticalProteins')
# Positive cases
.testEquals(getNCBIIdenticalProteins('AHA80958', format = 'ids'), list(c('WP_063864654.1','AHA80958.1', "EKD8974449.1", "EKD8979565.1")))
.testEquals(getNCBIIdenticalProteins(c('AHA80958', 'AHA80958'), format = 'ids'), list(c('WP_063864654.1','AHA80958.1', "EKD8974449.1", "EKD8979565.1"),
                                                                                     c('WP_063864654.1','AHA80958.1', "EKD8974449.1", "EKD8979565.1")))
Sys.sleep(3)
dummyDf <- data.frame(Id=c(45721358, 45721358, 45721358, 45721358),
                      Source=c('RefSeq','INSDC','INSDC','INSDC'),
                      Nucleotide.Accession=c('NG_050043.1','KF513177.1', 'ABJLVL010000001.1', 'ABJLVL010000113.1'),
                      Start=c(1,1, 124981, 1755),
                      Stop=c(861,861, 125841, 2615),Strand=c('+','+','-','+'),
                      Protein=c('WP_063864654.1','AHA80958.1', 'EKD8974449.1', 'EKD8979565.1'),
                      Protein.Name=c('class A beta-lactamase SHV-172','beta-lactamase SHV-172', 'class A beta-lactamase SHV-172',
                                     'class A beta-lactamase SHV-172'),
                      Organism=c('Klebsiella pneumoniae','Klebsiella pneumoniae','Klebsiella pneumoniae','Klebsiella pneumoniae'),
                      Strain=c(845332,845332, NA, NA),Assembly=c('', '','GCA_026265195.1','GCA_026265195.1'))
.testEquals(getNCBIIdenticalProteins('AHA80958', format = 'dataframe'), dummyDf)
.testEquals(getNCBIIdenticalProteins(c('AHA80958','AHA80958'), format = 'dataframe'), list(dummyDf, dummyDf))
# Incorrect format request case
Sys.sleep(3)
RUnit::checkException(getNCBIIdenticalProteins('AHA80958', format='test'))
# No identical proteins found case
# .testEquals(getNCBIIdenticalProteins('test'), NA)

### Test .getNCBI2UniProtBatch
# message('Testing .getNCBI2UniProtBatch')
# # Positive cases
# dummyDf <- getNCBIIdenticalProteins('CAA76794', format = 'dataframe')
# .testEquals(.getNCBI2UniProtBatch(dummyDf[which(dummyDf$Source=='INSDC'),] ,'EMBL-GenBank-DDBJ_CDS' ),
#             data.frame(From=c('AAC06040.1','AAN61404.1','CAA76794.1'),
#             Entry=c('O68642','Q8GP08','Q9R747')))
# # No translations cases
# dummyDf <- getNCBIIdenticalProteins('WP_041918279', format = 'dataframe')
# .testEquals(.getNCBI2UniProtBatch(dummyDf[which(dummyDf$Source=='RefSeq'),] ,'RefSeq_Protein' ),
#             data.frame(From=logical(),Entry=logical()))
# .testEquals(.getNCBI2UniProtBatch(dummyDf[which(dummyDf$Source=='INSDC'),] ,'EMBL-GenBank-DDBJ_CDS' ),
#             data.frame(From=logical(),Entry=logical()))

### Test .getNCBI2UniProtIP
# message('Testing .getNCBI2UniProtIP')
# # Positive case
# .testEquals(.getNCBI2UniProtIP('WP_010896559.1'), c('Q7AJZ0','Q9RC37'))
# # No identical proteins found case
# .testEquals(.getNCBI2UniProtIP('test'), character(0))
# # Identical proteins found, but no possible translation
# .testEquals(.getNCBI2UniProtIP('AMR06225.1'), character(0))

### Test getNCBIProtein2UniProt
message('Testing getNCBIProtein2UniProt')
# Positive cases
.testEquals(getNCBIProtein2UniProt('WP_010896559.1'), c('Q7AJZ0'))
.testEquals(getNCBIProtein2UniProt(c('WP_010896559.1','WP_010896559.1')), c('Q7AJZ0', 'Q7AJZ0'))
# .testEquals(getNCBIProtein2UniProt('WP_010896559.1', exhaustiveMapping = TRUE), c('Q7AJZ0','Q9RC37'))
# .testEquals(getNCBIProtein2UniProt('WP_010896559.1', exhaustiveMapping = TRUE, detailedMapping = TRUE), list('DT'=c('Q7AJZ0', 'Q9RC37'),'1.0'=c('Q7AJZ0', 'Q9RC37')))
# .testEquals(getNCBIProtein2UniProt('WP_010896559.1', exhaustiveMapping = FALSE, detailedMapping = TRUE), list('DT'=c('Q7AJZ0')))
# .testEquals(getNCBIProtein2UniProt('WP_039189232.1'), c('A0A2A2MC99'))
# .testEquals(getNCBIProtein2UniProt('WP_039189232.1', exhaustiveMapping = TRUE, byIdenticalProteins = FALSE), c('A0A2A2MC99','Q9K351'))
# .testEquals(getNCBIProtein2UniProt('WP_039189232.1', byIdenticalProteins = FALSE), c('A0A2A2MC99'))
# No translation case
.testEquals(getNCBIProtein2UniProt('WP_188331862.1'), NA)
.testEquals(getNCBIProtein2UniProt('WP_188331862.1', detailedMapping = TRUE), NA)
# NCBI Protein ID not registered case
# .testEquals(getNCBIProtein2UniProt('test'), NA)

### Test getNCBINucleotide2UniProt
message('Testing getNCBINucleotide2UniProt')
# Positive case
.testEquals(getNCBINucleotide2UniProt('AY536519'), c('Q6QJ79'))
.testEquals(getNCBINucleotide2UniProt(c('AY536519', 'AY536519')), c('Q6QJ79','Q6QJ79'))
.testEquals(getNCBINucleotide2UniProt('AY536519', detailedMapping = TRUE), list('1.0'=c('Q6QJ79')))
# .testEquals(getNCBINucleotide2UniProt('AY536519', exhaustiveMapping = TRUE), c('Q6QJ79','A0A7G1KXU2','A0A6I4WTI5','D0UY02'))
# .testEquals(getNCBINucleotide2UniProt('AY536519', exhaustiveMapping = TRUE, detailedMapping = TRUE), list('1.0'=c('Q6QJ79','A0A7G1KXU2','A0A6I4WTI5','D0UY02')))
# .testEquals(getNCBINucleotide2UniProt('AY536519', byIdenticalProteins = FALSE), character(0))
# No translation case
.testEquals(getNCBINucleotide2UniProt('Z21488'), NA)
.testEquals(getNCBINucleotide2UniProt('Z21488', detailedMapping = TRUE), NA)
# NCBI Nucleotide ID not registered case
# .testEquals(getNCBINucleotide2UniProt('test'), NA)

### Test getNCBIGene2UniProt
message('Testing getNCBIGene2UniProt')
# Positive case
.testEquals(getNCBIGene2UniProt('76524190'), c('A0A1W6DPG3'))
.testEquals(getNCBIGene2UniProt(c('76524190','76524190')), c('A0A1W6DPG3','A0A1W6DPG3'))
.testEquals(getNCBIGene2UniProt('76524190', detailedMapping = TRUE), list('1.0'=c('A0A1W6DPG3')))
.testEquals(getNCBIGene2UniProt('76524190', byIdenticalProteins =  FALSE), NA)
# No translation case
.testEquals(getNCBIGene2UniProt('69412715'), NA)
.testEquals(getNCBIGene2UniProt('69412715', detailedMapping = TRUE), NA)
# NCBI Gene ID not registered case
# .testEquals(getNCBIGene2UniProt('test'), NA)

##########################
# NCBI databases to KEGG #
##########################

### Test .getNCBI2KEGGDT
# message('Testing .getNCBI2KEGGDT')
# # Positive cases
# .testEquals(.getNCBI2KEGGDT('AAC44793.1'), c('ag:AAC44793'))
# .testEquals(.getNCBI2KEGGDT('AAO66446.1'), c('ag:AAO66446'))
# .testEquals(.getNCBI2KEGGDT('NP_059345'), c('hsa:10458'))
# .testEquals(.getNCBI2KEGGDT('AAG58814'), c('ece:Z5100'))
# # No translation case
# .testEquals(.getNCBI2KEGGDT('CAD69003.1'), character(0))

### Test .getNCBI2KEGGTUP
# message('Testing .getNCBI2KEGGTUP')
# # Positive cases
# # .testEquals(.getNCBI2KEGGTUP('WP_010896559.1', 'protein', exhaustiveMapping = TRUE, detailedMapping = TRUE, byIdenticalProteins = FALSE),
# #             list('DT'=c('bha:BH0380'), '0.5'=c("ag:AAA22599","vpn:A21D_00889","ag:AAP74657","ddh:Desde_1644","bcl:ABC3508","scib:HUG20_05815","bha:BH0380"),
# #                  '1.0'=c('bha:BH0380'), '0.9'=c('bha:BH0380')))
# # .testEquals(.getNCBI2KEGGTUP('WP_010896559.1', 'protein', exhaustiveMapping = TRUE, detailedMapping = FALSE, byIdenticalProteins = FALSE),
# #             c("bha:BH0380","ag:AAA22599","vpn:A21D_00889","ag:AAP74657","ddh:Desde_1644","bcl:ABC3508","scib:HUG20_05815"))
# .testEquals(.getNCBI2KEGGTUP('WP_010896559.1', 'protein', exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = FALSE),
#             c("bha:BH0380"))
# .testEquals(.getNCBI2KEGGTUP('WP_010896559.1', 'protein', exhaustiveMapping = FALSE, detailedMapping = TRUE, byIdenticalProteins = FALSE),
#             list('DT'=c("bha:BH0380")))
# # No translation case
# .testEquals(.getNCBI2KEGGTUP('WP_188331862.1', 'protein'), character(0))
# .testEquals(.getNCBI2KEGGTUP('WP_188331862.1', 'protein', detailedMapping = TRUE), list())


### Test getNCBIProtein2KEGG
message('Testing getNCBIProtein2KEGG')
# Positive cases
.testEquals(getNCBIProtein2KEGG('WP_010896559.1', exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = FALSE),
            c("bha:BH0380"))
.testEquals(getNCBIProtein2KEGG(c('WP_010896559.1','WP_010896559.1'), exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = FALSE),
           c("bha:BH0380", "bha:BH0380"))
.testEquals(getNCBIProtein2KEGG('WP_010896559.1', exhaustiveMapping = FALSE, detailedMapping = TRUE, byIdenticalProteins = FALSE),
            list('DT'=c("bha:BH0380")))
.testEquals(getNCBIProtein2KEGG(c('WP_010896559.1','WP_010896559.1'), exhaustiveMapping = FALSE, detailedMapping = TRUE, byIdenticalProteins = FALSE),
            c('DT'=c("bha:BH0380"),'DT'=c("bha:BH0380")))
# .testEquals(getNCBIProtein2KEGG('WP_010896559.1',exhaustiveMapping = TRUE, detailedMapping = TRUE, byIdenticalProteins = FALSE),
#             list('DT'=c('bha:BH0380'), '0.5'=c("ag:AAA22599","vpn:A21D_00889","ag:AAP74657","ddh:Desde_1644","bcl:ABC3508","scib:HUG20_05815","bha:BH0380"),
#                  '1.0'=c('bha:BH0380'), '0.9'=c('bha:BH0380')))
# .testEquals(getNCBIProtein2KEGG('WP_010896559.1',exhaustiveMapping = TRUE, detailedMapping = FALSE, byIdenticalProteins = FALSE),
#             c("bha:BH0380","ag:AAA22599","vpn:A21D_00889","ag:AAP74657","ddh:Desde_1644","bcl:ABC3508","scib:HUG20_05815"))
# No translation case
.testEquals(getNCBIProtein2KEGG('WP_188331862.1'), NA)
.testEquals(getNCBIProtein2KEGG('WP_188331862.1', detailedMapping = TRUE), NA)
# NCBI Protein ID not registered case
# .testEquals(getNCBIProtein2KEGG('test'), NA)

### Test getNCBINucleotide2KEGG
message('Testing getNCBINucleotide2KEGG')
# Positive cases
.testEquals(getNCBINucleotide2KEGG('AY536519', exhaustiveMapping = FALSE, detailedMapping = FALSE),
            c("ag:AAS48620"))
.testEquals(getNCBINucleotide2KEGG(c('AY536519','AY536519'), exhaustiveMapping = FALSE, detailedMapping = FALSE),
           c("ag:AAS48620", "ag:AAS48620"))
# .testEquals(getNCBINucleotide2KEGG('NZ_CP059690', exhaustiveMapping = FALSE, detailedMapping = TRUE, byIdenticalProteins = FALSE),
#             list('1.0'=c("apa:APP7_0365")))
# No translation case
.testEquals(getNCBINucleotide2KEGG('AY536519', byIdenticalProteins = FALSE), NA)
.testEquals(getNCBINucleotide2KEGG('AY536519', detailedMapping = TRUE, byIdenticalProteins = FALSE), NA)
# NCBI Nucleotide ID not registered case
# .testEquals(getNCBINucleotide2KEGG('test'), NA)

### Test getNCBIGene2KEGG
message('Testing getNCBIGene2KEGG')
# Positive cases
.testEquals(getNCBIGene2KEGG('76524190', exhaustiveMapping = FALSE, detailedMapping = FALSE),
            c("cet:B8281_09025"))
.testEquals(getNCBIGene2KEGG(c('76524190','76524190'), exhaustiveMapping = FALSE, detailedMapping = FALSE),
           c("cet:B8281_09025","cet:B8281_09025"))
# .testEquals(getNCBIGene2KEGG('76524190', exhaustiveMapping = FALSE, detailedMapping = TRUE),
#             list('DT'=c("cet:B8281_09025")))
# No translation case
.testEquals(getNCBIGene2KEGG('76524190', byIdenticalProteins = FALSE, bySimilarGenes = FALSE), NA)
.testEquals(getNCBIGene2KEGG('76524190', detailedMapping = TRUE, byIdenticalProteins = FALSE, bySimilarGenes = FALSE), NA)
# NCBI Gene ID not registered case
# .testEquals(getNCBIGene2KEGG('test'), NA)

##########################
# NCBI databases to CARD #
##########################

### Test getNCBIProtein2CARD
message('Testing getNCBIProtein2CARD')
# Positive cases
.testEquals(getNCBIProtein2CARD('AAK64581'), c('ARO:3004568'))
.testEquals(getNCBIProtein2CARD(c('AAK64581','AAK64581')), c('ARO:3004568', 'ARO:3004568'))
.testEquals(getNCBIProtein2CARD('AAK64581.1'), c('ARO:3004568'))
.testEquals(getNCBIProtein2CARD('CAC81324'), c('ARO:3003015'))
.testEquals(getNCBIProtein2CARD('CAC81324.1'), c('ARO:3003015'))
# No translation case
.testEquals(getNCBIProtein2CARD('EAT32321.1'), NA)
# NCBI ID not registered case
# .testEquals(getNCBIProtein2CARD('test'), NA)

### Test getNCBINucleotide2CARD
message('Testing getNCBINucleotide2CARD')
# Positive cases
.testEquals(getNCBINucleotide2CARD('AY034138'), c('ARO:3004568'))
.testEquals(getNCBINucleotide2CARD(c('AY034138','AY034138')), c('ARO:3004568','ARO:3004568'))
.testEquals(getNCBINucleotide2CARD('AY034138.1'), c('ARO:3004568'))
.testEquals(getNCBINucleotide2CARD('AJ310778'), c('ARO:3003015'))
.testEquals(getNCBINucleotide2CARD('AJ310778.1'), c('ARO:3003015'))
# No translation case
.testEquals(getNCBINucleotide2CARD('JH725437.1'), NA)
# NCBI ID not registered case
# .testEquals(getNCBINucleotide2CARD('test'), NA)

### Test getNCBIGene2CARD
message('Testing getNCBIGene2CARD')
# Positive case
.testEquals(getNCBIGene2CARD('3510143'), c('ARO:3003942'))
.testEquals(getNCBIGene2CARD(c('3510143','3510143')), c('ARO:3003942','ARO:3003942'))
# .testEquals(getNCBIGene2CARD('WP_001082319'), NA)
# NCBI ID not registered case
# .testEquals(getNCBIGene2CARD('test'), NA)


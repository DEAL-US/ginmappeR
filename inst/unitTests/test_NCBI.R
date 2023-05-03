library(RUnit)
library(ginmappeR)
library('UniProt.ws')
library('KEGGREST')
library('httr')
library('rentrez')
library('XML')
library('pkgfilecache')
# source('../R/NCBIFunctions.R')
# source('../R/UniProtFunctions.R')
# source('../R/utilsFunctions.R')

#####################################
# NCBI databases inter-translations #
#####################################

### Test getNCBIGene2NCBIProtein

# Positive case
checkEquals(getNCBIGene2NCBIProtein('76524190'), c('WP_001082319'))
checkEquals(getNCBIGene2NCBIProtein('76524190', TRUE), c('WP_001082319'))
# ID not registered case
checkException(getNCBIGene2NCBIProtein('test'))
# No translation case
checkEquals(getNCBIGene2NCBIProtein('WP_001082319'), character(0))

### Test getNCBIProtein2NCBIGene

# Positive case
checkEquals(getNCBIProtein2NCBIGene('CAA79696'), c('1272'))
checkEquals(getNCBIProtein2NCBIGene('CAA79696', TRUE), c('1272'))
# ID not registered case
checkException(getNCBIProtein2NCBIGene('test'))
# No translation case
checkEquals(getNCBIProtein2NCBIGene('WP_011997479'), character(0))

### Test getNCBIProtein2NCBINucleotide

# Positive case
checkEquals(getNCBIProtein2NCBINucleotide('CAA79696'), c('Z21488'))
checkEquals(getNCBIProtein2NCBINucleotide('CAA79696', TRUE), c('Z21488'))
# ID not registered case
checkException(getNCBIProtein2NCBINucleotide('test'))

### Test getNCBINucleotide2NCBIProtein

# # Positive case
checkEquals(getNCBINucleotide2NCBIProtein('AY536519'), c('AAS48620'))
checkEquals(getNCBINucleotide2NCBIProtein('AY536519', TRUE), c('AAS48620'))
# ID not registered case
checkException(getNCBINucleotide2NCBIProtein('test'))

### Test getNCBIGene2NCBINucleotide

# Positive case
checkEquals(getNCBIGene2NCBINucleotide('76524190'), c('NZ_CP059690'))
checkEquals(getNCBIGene2NCBINucleotide('76524190', TRUE), c('NZ_CP059690'))
# ID not registered case
checkException(getNCBIGene2NCBINucleotide('test'))
# No translation case
checkEquals(getNCBIGene2NCBINucleotide('WP_001082319'), character(0))

### Test getNCBINucleotide2NCBIGene

# Positive case
checkEquals(getNCBINucleotide2NCBIGene('Z21488'), c('1272'))
checkEquals(getNCBINucleotide2NCBIGene('Z21488', TRUE), c('1272'))
# ID not registered case
checkException(getNCBINucleotide2NCBIGene('test'))
# No translation case
checkEquals(getNCBINucleotide2NCBIGene('KF513177'), character(0))

#############################
# NCBI databases to UniProt #
#############################

### Test .getNCBI2UniProtDT

# Positive cases
checkEquals(.getNCBI2UniProtDT('AEJ08681'), c('F8TCS6'))
checkEquals(.getNCBI2UniProtDT('AEJ08681 AGQ48857.1 CAA79696'), c('F8TCS6','S5FUH0','Q12860'))
checkEquals(.getNCBI2UniProtDT('CAA79696'), c('Q12860'))
# ID not registered case
checkEquals(.getNCBI2UniProtDT('test'), character(0))

### Test getNCBIIdenticalProteins

# Positive cases
checkEquals(getNCBIIdenticalProteins('AHA80958', format = 'ids'), c('WP_063864654.1','AHA80958.1'))
dummyDf <- data.frame(Id=c(45721358, 45721358), Source=c('RefSeq','INSDC'),
                      Nucleotide.Accession=c('NG_050043.1','KF513177.1'),
                      Start=c(1,1),Stop=c(861,861),Strand=c('+','+'),
                      Protein=c('WP_063864654.1','AHA80958.1'),
                      Protein.Name=c('class A beta-lactamase SHV-172','beta-lactamase SHV-172'),
                      Organism=c('Klebsiella pneumoniae','Klebsiella pneumoniae'),
                      Strain=c(845332,845332),Assembly=c(NA,NA))
checkEquals(getNCBIIdenticalProteins('AHA80958', format = 'dataframe'), dummyDf)
# Incorrect format request case
checkException(getNCBIIdenticalProteins('AHA80958', format='test'))
# No identical proteins found case
checkEquals(getNCBIIdenticalProteins('test'), character(0))

### Test .getNCBI2UniProtBatch

# Positive cases
dummyDf <- getNCBIIdenticalProteins('CAA76794', format = 'dataframe')
checkEquals(.getNCBI2UniProtBatch(dummyDf[which(dummyDf$Source=='RefSeq'),] ,'RefSeq_Protein' ),
            data.frame(From=c('WP_063864899.1'), Entry=c('Q8GP08')))
checkEquals(.getNCBI2UniProtBatch(dummyDf[which(dummyDf$Source=='INSDC'),] ,'EMBL-GenBank-DDBJ_CDS' ),
            data.frame(From=c('AAC06040.1','AAN61404.1','CAA76794.1'),
            Entry=c('O68642','Q8GP08','Q9R747')))
# No translations cases
dummyDf <- getNCBIIdenticalProteins('WP_041918279', format = 'dataframe')
checkEquals(.getNCBI2UniProtBatch(dummyDf[which(dummyDf$Source=='RefSeq'),] ,'RefSeq_Protein' ),
            data.frame(From=logical(),Entry=logical()))
checkEquals(.getNCBI2UniProtBatch(dummyDf[which(dummyDf$Source=='INSDC'),] ,'EMBL-GenBank-DDBJ_CDS' ),
            data.frame(From=logical(),Entry=logical()))

### Test .getNCBI2UniProtIP

# Positive case
checkEquals(.getNCBI2UniProtIP('WP_010896559.1'), c('Q7AJZ0','Q9RC37'))
# No identical proteins found case
checkEquals(.getNCBI2UniProtIP('test'), character(0))
# Identical proteins found, but no possible translation
checkEquals(.getNCBI2UniProtIP('AMR06225.1'), character(0))


### Test getNCBIProtein2UniProt

# Positive cases
checkEquals(getNCBIProtein2UniProt('WP_010896559.1'), c('Q7AJZ0'))
checkEquals(getNCBIProtein2UniProt('WP_010896559.1', TRUE), c('Q7AJZ0','Q9RC37'))
checkEquals(getNCBIProtein2UniProt('WP_039189232.1'), c('A0A2A2MC99','Q9K351', 'A0A141RJY1'))
checkEquals(getNCBIProtein2UniProt('WP_039189232.1', TRUE, byIdenticalProteins = FALSE), c('A0A2A2MC99','Q9K351'))
checkEquals(getNCBIProtein2UniProt('WP_039189232.1', byIdenticalProteins = FALSE), c('A0A2A2MC99'))
# No translation case
checkEquals(getNCBIProtein2UniProt('WP_188331862.1'), character(0))
# NCBI Protein ID not registered case
checkException(getNCBIProtein2UniProt('test'))

### Test getNCBINucleotide2UniProt

# Positive case
checkEquals(getNCBINucleotide2UniProt('AY536519'), c('Q6QJ79'))
checkEquals(getNCBINucleotide2UniProt('AY536519', TRUE), c('Q6QJ79','A0A7G1KXU2','A0A6I4WTI5','D0UY02'))
checkEquals(getNCBINucleotide2UniProt('AY536519', byIdenticalProteins = FALSE), character(0))
# No translation case
checkEquals(getNCBINucleotide2UniProt('Z21488'), character(0))
# NCBI Nucleotide ID not registered case
checkException(getNCBINucleotide2UniProt('test'))

### Test getNCBIGene2UniProt

# Positive case
checkEquals(getNCBIGene2UniProt('76524190'), c('A0A1W6DPG3'))
checkEquals(getNCBIGene2UniProt('76524190', byIdenticalProteins =  FALSE), character(0))
# No translation case
checkEquals(getNCBIGene2UniProt('69412715'), character(0))
# NCBI Gene ID not registered case
checkException(getNCBIGene2UniProt('test'))

##########################
# NCBI databases to KEGG #
##########################

### Test .getNCBI2KEGGDT

# Positive cases
checkEquals(.getNCBI2KEGGDT('AAC44793.1'), c('ag:AAC44793'))
checkEquals(.getNCBI2KEGGDT('AAO66446.1'), c('ag:AAO66446'))
checkEquals(.getNCBI2KEGGDT('NP_059345'), c('hsa:10458'))
checkEquals(.getNCBI2KEGGDT('AAG58814'), c('ece:Z5100'))
# No translation case
checkEquals(.getNCBI2KEGGDT('CAD69003.1'), character(0))

### Test .getNCBI2KEGGTUP

# Positive cases
checkEquals(.getNCBI2KEGGTUP('WP_010896559.1', 'protein', exhaustiveMapping = TRUE, detailedMapping = TRUE, byIdenticalProteins = FALSE),
            list('DT'=c('bha:BH0380'), '0.5'=c("ag:AAA22599","vpn:A21D_00889","ag:AAP74657","ddh:Desde_1644","bcl:ABC3508","scib:HUG20_05815","bha:BH0380"),
                 '1.0'=c('bha:BH0380'), '0.9'=c('bha:BH0380')))
checkEquals(.getNCBI2KEGGTUP('WP_010896559.1', 'protein', exhaustiveMapping = TRUE, detailedMapping = FALSE, byIdenticalProteins = FALSE),
            c("bha:BH0380","ag:AAA22599","vpn:A21D_00889","ag:AAP74657","ddh:Desde_1644","bcl:ABC3508","scib:HUG20_05815"))
checkEquals(.getNCBI2KEGGTUP('WP_010896559.1', 'protein', exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = FALSE),
            c("bha:BH0380"))
checkEquals(.getNCBI2KEGGTUP('WP_010896559.1', 'protein', exhaustiveMapping = FALSE, detailedMapping = TRUE, byIdenticalProteins = FALSE),
            list('DT'=c("bha:BH0380")))
# No translation case
checkEquals(.getNCBI2KEGGTUP('WP_188331862.1', 'protein'), character(0))
checkEquals(.getNCBI2KEGGTUP('WP_188331862.1', 'protein', detailedMapping = TRUE), list())


### Test getNCBIProtein2KEGG

# Positive cases
checkEquals(getNCBIProtein2KEGG('WP_010896559.1', exhaustiveMapping = FALSE, detailedMapping = FALSE, byIdenticalProteins = FALSE),
            c("bha:BH0380"))
checkEquals(getNCBIProtein2KEGG('WP_010896559.1', exhaustiveMapping = FALSE, detailedMapping = TRUE, byIdenticalProteins = FALSE),
            list('DT'=c("bha:BH0380")))
checkEquals(getNCBIProtein2KEGG('WP_010896559.1',exhaustiveMapping = TRUE, detailedMapping = TRUE, byIdenticalProteins = FALSE),
            list('DT'=c('bha:BH0380'), '0.5'=c("ag:AAA22599","vpn:A21D_00889","ag:AAP74657","ddh:Desde_1644","bcl:ABC3508","scib:HUG20_05815","bha:BH0380"),
                 '1.0'=c('bha:BH0380'), '0.9'=c('bha:BH0380')))
checkEquals(getNCBIProtein2KEGG('WP_010896559.1',exhaustiveMapping = TRUE, detailedMapping = FALSE, byIdenticalProteins = FALSE),
            c("bha:BH0380","ag:AAA22599","vpn:A21D_00889","ag:AAP74657","ddh:Desde_1644","bcl:ABC3508","scib:HUG20_05815"))
# No translation case
checkEquals(getNCBIProtein2KEGG('WP_188331862.1'), character(0))
checkEquals(getNCBIProtein2KEGG('WP_188331862.1', detailedMapping = TRUE), list())

### Test getNCBINucleotide2KEGG

# Positive cases

checkEquals(getNCBINucleotide2KEGG('AY536519', exhaustiveMapping = FALSE, detailedMapping = FALSE),
            c("ag:AAS48620"))
checkEquals(getNCBINucleotide2KEGG('NZ_CP059690', exhaustiveMapping = FALSE, detailedMapping = TRUE, byIdenticalProteins = FALSE),
            list('1.0'=c("apa:APP7_0365","apa:APP7_1820")))
# No translation case
checkEquals(getNCBINucleotide2KEGG('AY536519', byIdenticalProteins = FALSE), character(0))
checkEquals(getNCBINucleotide2KEGG('AY536519', detailedMapping = TRUE, byIdenticalProteins = FALSE), list())

### Test getNCBIGene2KEGG

# Positive cases
checkEquals(getNCBIGene2KEGG('76524190', exhaustiveMapping = FALSE, detailedMapping = FALSE),
            c("cet:B8281_09025"))
checkEquals(getNCBIGene2KEGG('76524190', exhaustiveMapping = FALSE, detailedMapping = TRUE),
            list('DT'=c("cet:B8281_09025")))
# No translation case
checkEquals(getNCBIGene2KEGG('76524190', byIdenticalProteins = FALSE, bySimilarGenes = FALSE), character(0))
checkEquals(getNCBIGene2KEGG('76524190', detailedMapping = TRUE, byIdenticalProteins = FALSE, bySimilarGenes = FALSE), list())

##########################
# NCBI databases to CARD #
##########################

### Test getNCBIProtein2CARD

# Positive cases
checkEquals(getNCBIProtein2CARD('AAK64581'), c('ARO:3004568'))
checkEquals(getNCBIProtein2CARD('AAK64581.1'), c('ARO:3004568'))
checkEquals(getNCBIProtein2CARD('CAC81324'), c('ARO:3003015'))
checkEquals(getNCBIProtein2CARD('CAC81324.1'), c('ARO:3003015'))
# No translation case
checkEquals(getNCBIProtein2CARD('EAT32321.1'), character(0))

### Test getNCBINucleotide2CARD

# Positive cases
checkEquals(getNCBINucleotide2CARD('AY034138'), c('ARO:3004568'))
checkEquals(getNCBINucleotide2CARD('AY034138.1'), c('ARO:3004568'))
checkEquals(getNCBINucleotide2CARD('AJ310778'), c('ARO:3003015'))
checkEquals(getNCBINucleotide2CARD('AJ310778.1'), c('ARO:3003015'))
# No translation case
checkEquals(getNCBINucleotide2CARD('XM_001230804.1'), character(0))

### Test getNCBIGene2CARD

# Positive case
checkEquals(getNCBIGene2CARD('3510143'), c('ARO:3003942'))
checkEquals(getNCBIGene2CARD('WP_001082319'), character(0))













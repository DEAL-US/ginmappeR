library(RUnit)
source('../R/NCBIFunctions.R')
source('../R/utilsFunctions.R')
library('UniProt.ws')
library('rentrez')
library('KEGGREST')
library('XML')

#####################################
# NCBI databases inter-translations #
#####################################

### Test getNCBIGene2NCBIProtein

# Positive case
checkEquals(getNCBIGene2NCBIProtein('76524190'), list('WP_001082319'))
# ID not registered case
checkException(getNCBIGene2NCBIProtein('test'))
# No translation case
checkException(getNCBIGene2NCBIProtein('WP_001082319'))

### Test getNCBIProtein2NCBIGene

# Positive case
checkEquals(getNCBIProtein2NCBIGene('CAA79696'), list('1272'))
# ID not registered case
checkException(getNCBIProtein2NCBIGene('test'))
# No translation case
checkException(getNCBIProtein2NCBIGene('WP_011997479'))

### Test getNCBIProtein2NCBINucleotide

# Positive case
checkEquals(getNCBIProtein2NCBINucleotide('CAA79696'), list('Z21488'))
# ID not registered case
checkException(getNCBIProtein2NCBINucleotide('test'))

### Test getNCBINucleotide2NCBIProtein

# # Positive case
checkEquals(getNCBINucleotide2NCBIProtein('AY536519'), list('AAS48620'))
# ID not registered case
checkException(getNCBINucleotide2NCBIProtein('test'))

### Test getNCBIGene2NCBINucleotide

# Positive case
checkEquals(getNCBIGene2NCBINucleotide('76524190'), list('NZ_CP059690'))
# ID not registered case
checkException(getNCBIGene2NCBINucleotide('test'))
# No translation case
checkException(getNCBIGene2NCBINucleotide('WP_001082319'))

### Test getNCBINucleotide2NCBIGene

# Positive case
checkEquals(getNCBINucleotide2NCBIGene('Z21488'), list('1272'))
# ID not registered case
checkException(getNCBINucleotide2NCBIGene('test'))
# No translation case
checkException(getNCBINucleotide2NCBIGene('KF513177'))

#####################################
# NCBI databases auxiliar functions #
#####################################

### Test .getNCBI2UniProtDT

# Positive case
checkEquals(.getNCBI2UniProtDT('AEJ08681'), list('F8TCS6'))
checkEquals(.getNCBI2UniProtDT('AEJ08681 AGQ48857.1 CAA79696'), list('F8TCS6','S5FUH0','Q12860'))
checkEquals(.getNCBI2UniProtDT('CAA79696'), list('Q12860'))
# ID not registered case
checkEquals(.getNCBI2UniProtDT('test'), list())

### Test getNCBIIdenticalProteins

# Positive cases
checkEquals(getNCBIIdenticalProteins('AHA80958', format = 'ids'), list('WP_063864654.1','AHA80958.1'))
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
checkException(getNCBIIdenticalProteins('test'))

### Test .getNCBI2UniProtBatch

# Positive case
dummyDf <- getNCBIIdenticalProteins('CAA76794', format = 'dataframe')
checkEquals(.getNCBI2UniProtBatch(dummyDf[which(dummyDf$Source=='RefSeq'),] ,'RefSeq_Protein' ),
            data.frame(From=c('WP_063864899.1'), Entry=c('Q8GP08')))
checkEquals(.getNCBI2UniProtBatch(dummyDf[which(dummyDf$Source=='INSDC'),] ,'EMBL-GenBank-DDBJ_CDS' ),
            data.frame(From=c('AAC06040.1','AAN61404.1','CAA76794.1'),
            Entry=c('O68642','Q8GP08','Q9R747')))
# No translations case
dummyDf <- getNCBIIdenticalProteins('WP_041918279', format = 'dataframe')
checkEquals(.getNCBI2UniProtBatch(dummyDf[which(dummyDf$Source=='RefSeq'),] ,'RefSeq_Protein' ),
            data.frame(From=logical(),Entry=logical()))
checkEquals(.getNCBI2UniProtBatch(dummyDf[which(dummyDf$Source=='INSDC'),] ,'EMBL-GenBank-DDBJ_CDS' ),
            data.frame(From=logical(),Entry=logical()))

### Test .getNCBI2UniProtIP

# Positive case
checkEquals(.getNCBI2UniProtIP('WP_010896559.1'), list('Q7AJZ0','Q9RC37'))
# No identical proteins found case
checkEquals(.getNCBI2UniProtIP('test'), list())
# Identical proteins found, but no possible translation
checkEquals(.getNCBI2UniProtIP('AMR06225.1'), list())


### Test getNCBIProtein2UniProt

# Positive cases
checkEquals(getNCBIProtein2UniProt('WP_010896559.1'), list('Q7AJZ0','Q9RC37'))
checkEquals(getNCBIProtein2UniProt('WP_039189232.1'), list('A0A2A2MC99','Q9K351', 'A0A141RJY1'))
# No translation case
checkException(getNCBIProtein2UniProt('WP_188331862.1'))
# NCBI Protein ID not registered case
checkException(getNCBIProtein2UniProt('test'))

### Test getNCBINucleotide2UniProt

# Positive cases
checkEquals(getNCBINucleotide2UniProt('AY536519'), list('Q6QJ79','A0A7G1KXU2','A0A6I4WTI5','D0UY02'))
# No translation case
checkException(getNCBINucleotide2UniProt('Z21488'))
# NCBI Nucleotide ID not registered case
checkException(getNCBINucleotide2UniProt('test'))

### Test getNCBIGene2UniProt

# Positive cases
checkEquals(getNCBIGene2UniProt('76524190'), list('A0A1W6DPG3'))
# No translation case
checkException(getNCBIGene2UniProt('69412715'))
# NCBI Gene ID not registered case
checkException(getNCBIGene2UniProt('test'))


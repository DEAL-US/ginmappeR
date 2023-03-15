library(RUnit)

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








# # Code grabbed from https://www.rdocumentation.org/packages/RUnit/versions/0.4.32
#
# pkgname <- 'biomappeR'
# require(pkgname, quietly=TRUE, character.only=TRUE) || stop("package '", pkgname, "' not found")
#
# # How to determine which files to load (have to start with test_ and end with .R)
# # pattern <- "^test_.*\\.R$"
# pattern <- '.*'
#
# # Which functions to run. Have to start with 'test.'
# testFunctionRegexp = '.*'
#
# # Path to the unit tests folder in the package
# dir <- system.file('inst/tests', package=pkgname)
#
# # Define RUnit test suite
# suite <- RUnit:::defineTestSuite(name=paste(pkgname, "RUnit Tests"),
#                          dirs=dir,
#                          testFileRegexp=pattern,
#                          testFuncRegexp = testFunctionRegexp,
#                          rngKind="default",
#                          rngNormalKind="default")
#
# # Run tests
# result <- RUnit:::runTestSuite(suite)
#
# # Display result tests on the console
# printTextProtocol(result)

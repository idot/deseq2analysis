## for manual usage from pachage root while developing:
## must be pasted int console after test()
filesPath <- function(f){ paste("inst/testFiles/", f, sep="") }

# for tests:
filesPath <- function(f){ system.file("testFiles", f, package="deseq2analysis") }

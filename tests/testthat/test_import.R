context("import")


source("helper.R")


test_that("we can import an ensembl metadata file", {
  ensembl <- readEnsembl(filesPath("test.ensembl.GRCh38.genes.tab"))
  expect_equal(98, nrow(ensembl))
  expect_equal(8, ncol(ensembl))
  expect_equal("geneid", colnames(ensembl)[1])
})


test_that("we can import a counts file from featureCounts", {
   counts <- readCounts(filesPath("test.counts.tab"))
   expect_equal(98, nrow(counts))
   expect_equal(15, ncol(counts))
   expect_equal(as.character(83653:83667), colnames(counts))
})

test_that("we can import a counts file from featureCounts and it filters some genes", {
   counts <- readCounts(filesPath("test.counts.tab"), remove_genes=c("ENSG00000268020", "ENSG00000269732", "FOO"))
   expect_equal(96, nrow(counts))
   expect_equal(15, ncol(counts))
   expect_equal(as.character(83653:83667), colnames(counts))
})

test_that("we can get a metadata matrix from a grouping table", {
    grouping <- readGrouping(filesPath("test.grouping.tab"))
    expect_equal(15, nrow(grouping))
    expect_equal(2, ncol(grouping))
    expect_equal(as.character(83653:83667), grouping[,1])
})

test_that("we can ensure that the column names and the sample ids in the grouping are identical", {
  correct <- readGrouping(filesPath("test.grouping.tab"))
  wrong <- readGrouping(filesPath("test.grouping.wrong.tab"))
  counts <- readCounts(filesPath("test.counts.tab"))
  expect_true(checkCountColumnsDesignSamples(counts, correct))
  expect_error(checkCountColumnsDesignSamples(counts, wrong))
})


test_that("we can not get a metadata matrix from a grouping table and a sampleid vector that is sorted differently", {
  expect_error(readGrouping(filesPath("test.grouping.wrong.tab"),as.character(83653:83667)))
})

test_that("a deseq2 dataset can be created from a counts matrix and a grouping", {
     grouping <- readGrouping(filesPath("test.grouping.tab"))
     counts <- readCounts(filesPath("test.counts.tab"))
     dds <- deseqDataFromMatrix(counts, grouping)
     testthat::expect_is(dds, "DESeqDataSet")
     expect_equal(as.character(83653:83667), colnames(dds))
})

test_that("an error will be thrown when there is a difference between the counts columns and the grouping", {
  grouping <- readGrouping(filesPath("test.grouping.wrong.tab"))
  counts <- readCounts(filesPath("test.counts.tab"))
  expect_error(deseqDataFromMatrix(counts, grouping))
})

test_that("a deseq2 dataset with less columns can be created from a counts matrix and a grouping with filtered rows", {
  groupingAll <- readGrouping(filesPath("test.grouping.tab"))
  groupingFiltered <- groupingAll[groupingAll$group != "group2",]
  counts <- readCounts(filesPath("test.counts.tab"))
  dds <- deseqDataFromMatrix(counts, groupingFiltered)
  testthat::expect_is(dds, "DESeqDataSet")
  expect_equal(15, nrow(groupingAll))
  expect_gt( nrow(groupingAll), nrow(groupingFiltered))
  expect_equal(nrow(groupingFiltered), length(colnames(dds)))
})

test_that("we can get all the comparisons from a grouping as a two column tibble", {
  grouping <- readGrouping(filesPath("test.grouping.tab"))
  comparisons <- getComparisons(grouping)
  testthat::expect_is(comparisons, "tbl_df")
  expect_equal(10, nrow(comparisons))
  expect_equal(2, ncol(comparisons))
})

### to create test data
createTestData <- function(){
  filesPath <- function(f){ paste("inst/testFiles/", f, sep="") }
  grouping <- readGrouping(filesPath("test.grouping.tab"))
  counts <- readCounts(filesPath("test.counts.tab"))
  dds <- deseqDataFromMatrix(counts, grouping)
  saveRDS(dds, filesPath("test.dds.RDS"))
  dds.r <- DESeq2::DESeq(dds, betaPrior = TRUE)
  saveRDS(dds.r,filesPath("test.dds.r.RDS"))
  ### this is a special method that adds a comparison to metadata(deseq.result)
  deseq.result.m <- deseqResult(dds.r, "group2", "wt")
  saveRDS(deseq.result, filesPath("test.deseq.m.result.RDS"))
  deseq.result <- DESeq2::results(dds.r, contrast=c("group", "group2", "wt"))
  saveRDS(deseq.result, filesPath("test.deseq.result.RDS"))
}
####





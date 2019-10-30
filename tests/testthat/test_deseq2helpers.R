context("deseq2 helpers")

source("helper.R")

test_that("we can extract different names from a deseq result", {
  deseq.result <- readRDS(filesPath("test.deseq.result.RDS"))
  expect_equal("log2 fold change (MAP): group group2 vs wt", getComparisonFoldChange(deseq.result))
  expect_equal("group group2 vs wt", getComparisonTitle(deseq.result))
  expect_equal("group2_vs_wt", getComparisonString(deseq.result))
  expect_equal("group2 vs wt", getComparison(deseq.result))
})

test_that("we can the independent filtering threshold from the deseq.result", {
  deseq.result <- readRDS(filesPath("test.deseq.result.RDS"))
  expect_equivalent(0.0354, getFilterThreshold(deseq.result), tol = 1e-3 )
})

test_that("we get the normalised counts as a tibble", {
   co <- getCounts(readRDS(filesPath("test.dds.r.RDS")))
   expect_is(co, "tbl_df")
   expect_equal(98, length(co$geneid))
   expect_equal(16, ncol(co))
})

test_that("we get a deseq.result, but with the contrasts stored in metadata(comparison)", {
    dds.r <- readRDS(filesPath("test.dds.r.RDS"))
    dr <- deseqResult(dds.r, "group2", "group3")
    expect_is(dr,  'DESeqResults')
    expect_equal( c("group", "group2", "group3"), S4Vectors::metadata(dr)$comparison)
})

test_that("we get the normalised counts as a tibble but only from 2 groups", {
  dds.r <- readRDS(filesPath("test.dds.r.RDS"))
  grouping <- readGrouping(filesPath("test.grouping.tab"))
  co <- getCounts(dds.r, groups=c("group1","group2"), grouping=grouping)
  expect_is(co, "tbl_df")
  expect_equal(98, length(co$geneid))
  expect_equal(7, ncol(co))
  expect_equal(colnames(co), c(as.character(83656:83661),"geneid"))
})


test_that("we get the normalised counts as a tibble for all groups", {
  dds.r <- readRDS(filesPath("test.dds.r.RDS"))
  co <- getCounts(dds.r)
  expect_is(co, "tbl_df")
  expect_equal(98, length(co$geneid))
  expect_equal(16, ncol(co))
  expect_equal(colnames(co)[ncol(co)], "geneid")
})



test_that("we create a sorted tibble for saving and for plotting of the data", {
  ensembl <- readEnsembl(filesPath("test.ensembl.GRCh38.genes.tab"))
  deseq.result <- readRDS(filesPath("test.deseq.result.RDS"))
  comp <- toSortedTibble(deseq.result, ensembl, 0.01, ensembl_url)
  expect_equal(nrow(comp), nrow(deseq.result) )
  expect_equal(colnames(comp)[1],"geneid")
})

test_that("we create a sorted tibble for saving and for plotting of the data including normalised expression counts", {
  ensembl <- readEnsembl(filesPath("test.ensembl.GRCh38.genes.tab"))
  deseq.result <- readRDS(filesPath("test.deseq.m.result.RDS"))
  grouping <- readGrouping(filesPath("test.grouping.tab"))
  groups <- S4Vectors::metadata(deseq.result)$comparison[2:3]
  countsdata <- getCounts(readRDS(filesPath("test.dds.r.RDS")),groups)
  comp <- toSortedTibble(deseq.result, ensembl, 0.01, ensembl_url, countsdata)
  expect_equal(nrow(comp), nrow(deseq.result))
  #TODO: test width of comp
})




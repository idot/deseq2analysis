
#' reads counts table
#'
#' @param countpath path to aggregated counts file (column order is important!!!!) like from featureCounts
#' @param grouping grouping table
#' @param header2id converts column headers to ids, default: vbcf_bamname2id
#' @param removegenes vector of geneids to remove
#'
#' @return DESeq2 data object
#'
#' @export
readCountsToDeseq2 <- function(countpath, grouping, header2id=vbcf_bamname2id, remove_genes=NULL, metacols=6, skip=1){
  counts <- readCounts(countpath, header2id, remove_genes, metacols, skip)
  deseqDataFromMatrix(counts, grouping)
}


#' vbcf bam name to id
#'
#' extracts 5 consecutive digits from column name
#' @param column names to convert the ids from
#' @export
vbcf_bamname2id <- function(cols){
  stringr::str_extract(pattern="\\d{5}", string=cols)
}

#' 5 digits as id
#'
#' @param column names to convert the ids from
#' @export
digit5id <- function(cols){
  stringr::str_extract(pattern="\\d{5}", string=cols)
}

#' keep column header as is
#'
#' @param column names to convert the ids from
#' @export
keep <- function(cols){
  cols
}

#' reads ensembl info
#'
#' annotation should be a tab file with geneid, %gc, length columns
#' other columns are optional
#' creation of annotation with ensemblFA2tab.py
#' ensemblFA2tab input should consist of:
#'
#' zcat build.ncRNA.fa.gz build.cDNA.fa.gz > rna.fa  from ensembl
#' Ids and build should be identical to build.GTF
#'
#' @export
readEnsembl <- function(ensemblpath){
  ensgenes <- readr::read_tsv(ensemblpath) %>%
     dplyr::mutate( length_bin = cut(length, breaks=c(0,300,600,900,1500,3000,10000,1e9)), gc_bin = cut(gc, breaks=c(0,35,40,45,50,55,60,1)))
  ensgenes
}


#' reads count table
#'
#' @param countpath path to counts file
#' @param header2id function to convert column names to ids
#' @param remove_genes vector of ids to remove
#' @param metacols number of columns to remove because of metadata
#' @param skip number of rows to skip because of metadata#'
#'
#' @return matrix with remove_genes removed, columns renamed by header2id
#'
#' @export
readCounts <- function(countpath, header2id=vbcf_bamname2id, remove_genes=NULL, metacols=6, skip=1){
  counts <- readr::read_tsv(countpath,comment="",skip=skip)
  countsf <- counts %>% dplyr::filter(! .[[1]] %in% remove_genes)
  geneids <- countsf %>% .[[1]]
  countsfc <- countsf[,-c(1:metacols)]
  colnames(countsfc) <- header2id(colnames(countsfc))
  countsm <- as.matrix(countsfc)
  rownames(countsm) <- geneids
  countsm
}


#### TODO: base factor level???
#' gets grouping file to metadata matrix
#'
#' @param groupingpath path to grouping file with headers 1. column sampleId, 2. column group
#'
#' @export
readGrouping <- function(groupingpath){
  orig <- readr::read_tsv(groupingpath) %>% dplyr::mutate(sampleId=as.character(sampleId))
  grouping <- orig %>%
          dplyr::mutate(sampleId=as.character(sampleId)) %>%
          dplyr::group_by(group) %>%
          dplyr::arrange(sampleId) %>%
          dplyr::mutate(internal_replicate=1:dplyr::n()) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(group=factor(group))
  grouporder <- match(orig$sampleId, grouping$sampleId)
  groupingsorted <- grouping[grouporder,]
  as.data.frame(groupingsorted)
}



#' creates deseq2 DataSet with design by ~group from grouping table
#' countsMatrix columns can be filtered based on first column of grouping
#'
#' @export
deseqDataFromMatrix <- function(countsMatrix, grouping){
  countsMatrixS <- countsMatrix[,colnames(countsMatrix) %in% grouping[,1]]
  checkCountColumnsDesignSamples(countsMatrixS, grouping)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countsMatrixS,
                                colData = grouping,
                                design = ~ group)
  dds
}

#' checks that count matrix and grouing table have same sample order
#' grouping table must have sampleid in first column
#'
#' @export
checkCountColumnsDesignSamples <- function(countsMatrix, grouping){
  if(! all(colnames(countsMatrix) == grouping[,1])){
    err <- paste("problem with column order in count data and metadata\n", "   [" , paste(metadf[,1],collapse=", "), "]\n",  "   [", paste(countcolumns, collapse=","), "]\n", collapse=" ")
    stop(err)
  }
  TRUE
}

#' get the possible comparison combinations from the grouping table
#' @param grouping the grouping (with group as factor)
#'
#' @return a tibble with cond1 and cond2 as columns
#' @export
getComparisons <- function(grouping){
   allcomparisons <- data.frame(t(combn(unique(as.character(grouping$group)),2)),stringsAsFactors=FALSE)
   colnames(allcomparisons) <- c("cond1","cond2")
   tibble::as_tibble(allcomparisons)
}


#' read the comparisons table
#'
#' @param filepath path to comparisons table
#'
#' @export
readComparisonsTable <- function(filepath){
  readr::read_tsv(filepath)
}



#' gets normalised counts as a tibble from DeseqDataset
#'
#'
#' @return  a normalised counts tibble with geneid as last column
#'
#' @export
getCounts <- function(dds.r, groups=NULL, grouping=NULL){
  co <- DESeq2::counts(dds.r, normalized=TRUE)
  if(! is.null(groups) & ! is.null(grouping)){
     checkCountColumnsDesignSamples(dds.r,grouping)
     gi <-  which(grouping$group %in% groups)
     co <- co[,gi]
  }
  cot <- tibble::as_tibble(co)
  cot$geneid <- rownames(co)
  cot
}



#' deseq result for condition (default: group)
#'
#' group condition1 vs condition2
#'
#' normal deseq result with "contrast", adds condition names to metadata
#' could not find it, and parsing with getComparisonXXXX works but is brittle
#'
#' @param dds.r normalised DESeq dataset
#' @param condition1 the 1. condition
#' @param condition2 the 2. condition
#' @param condition (group to compare)
#'
#' @return the deseq2 result
#'
#' @export
deseqResult <- function(dds.r, condition1, condition2, condition="group"){
   comparison <- c(condition, condition1, condition2)
   deseq.result <- DESeq2::results(dds.r, contrast=comparison)
   S4Vectors::metadata(deseq.result)$comparison <- comparison
   deseq.result
}

#' converts deseq result to a sorted tibble joined with annotation
#'
#' annotation should be a tab file with geneid, %gc, length columns
#' other columns are optional
#' creation of annotation with ensemblFA2tab.py
#' ensemblFA2tab input should consist of:
#'
#' zcat build.ncRNA.fa.gz build.cDNA.fa.gz > rna.fa  from ensembl
#' Ids and build should be identical to build.GTF
#'
#' @param deseq.result deseq result dataset
#' @param ensembl annotation
#' @param countsdata the normalised expression values counts(dds.r, normalized=TRUE)
#' @param filterThreshold the mean value cutoff from deseq2 analysis
#'
#' @export
toSortedTibble <- function(deseq.result, ensembl, filterThreshold, countsdata=NULL){
    comp.all <-  as.data.frame(deseq.result) %>%
    dplyr::mutate(geneid = rownames(deseq.result), absFC = abs(log2FoldChange), meanFilter=baseMean < filterThreshold, url=ensembl_url(geneid), link=to_link(url, geneid))
    if(! is.null(countsdata)){
      comp.all <- comp.all %>% dplyr::left_join(countsdata)
    }
    comp.all <- comp.all %>%
    dplyr::left_join(ensembl, by=c(geneid="geneid")) %>%
    dplyr::arrange(padj, absFC) %>%
    dplyr::select(-absFC) %>%
    dplyr::select(geneid, dplyr::everything())
    if(nrow(deseq.result) != nrow(comp.all)){
      error <- paste("we lost some rows in toSortedTibble, should be left_join!: ", nrow(deseq.result), " -> ", nrow(comp.all))
      stop(error)
    }
    comp.all
}

#' extracts the independent filtering threshold from deseq result
#'
#' @export
getFilterThreshold <- function(deseq.result){
  S4Vectors::metadata(deseq.result)$filterThreshold
}

#' extracts the complete comparison string including log2FoldChange
#'
#' @export
getComparisonFoldChange <- function(deseq.result){
  meta <- S4Vectors::elementMetadata(deseq.result)
  lfcrow <- rownames(meta) == "log2FoldChange"
  meta$description[lfcrow]
}

#' extracts the comparison string
#' the comparison must be done on the group column
#'
#' @export
getComparison <- function(deseq.result){
  cfc <- getComparisonFoldChange(deseq.result)
  gsub(".*group (\\.*)", "\\1", cfc)
}

#' extracts the comparison
#'
#' @export
getComparisonTitle <- function(deseq.result){
  cn <- getComparisonFoldChange(deseq.result)
  trimws(strsplit(cn, ":")[[1]][2])
}

#' extracts a string from the comparison that is usable as a file path
#'
#' @export
getComparisonString <- function(deseq.result){
    # \\& \\ not working?
   #gsub("/| |\\||\\+|-|\\|\\&", "_", "as|d\\+f/sfd asdf& aer$|% ")
   REPLACER <- "/| |\\||\\+|-|\\|\\&"
   tit <- getComparisonTitle(deseq.result)
   gtit <- gsub(".*group (\\.*)", "\\1", tit)
   titu <- gsub(REPLACER, "_", gtit)
   titu
}


#' a meaningful explanation for the direction of the comparison
#' when log2 fold change > 0
#'
#' @export
getPosExplanation <- function(comparison){
   conds <- sapply(strsplit(comparison, "vs")[[1]], trimws)
   paste("log2 fold change > 0: ", "higher in ", conds[1], "lower in ", conds[2])
}

#' a meaningful explanation for the direction of the comparison
#' when log2 fold change < 0
#'
#' @export
getNegExplanation <- function(comparison){
  conds <- sapply(strsplit(comparison, "vs")[[1]], trimws)
  paste("log2 fold change < 0: ", "higher in ", conds[2], "lower in ", conds[1])
}

#' saving the result table in tab and xls
#'
#' @param comp.all a tibble with the data
#' @param pathbase the file name to save as will be appended with .diff.norm.[tab|xlsx]
#'
#' @export
saveTables <- function(comp.all, pathbase){
  readr::write_tsv(comp.all, paste(pathbase, ".diff.norm.tab",sep=""))
  writexl::write_xlsx(comp.all, paste(pathbase, ".diff.norm.xlsx",sep=""))
}

#' create an enseble url
#' @param create a url to ensembl gene id
#'
#' @export
ensembl_url <- function(ensemblid){
  sprintf("http://www.ensembl.org/id/%s",ensemblid)
}

#' create a link
#' @param url http://...
#' @param txt link text
#'
#' @export
to_link <- function(url, txt){
  paste('<a href="', url, '">', txt, '</a>', sep="")
}

#' extract row indices of significanlty deregulated genes
#' @param result list of deseq2 results
#' @param pcut adjusted p-value cutoff
#' @param lfc log fold change cutoff
#'
#' @export
extractSignif <- function(result, pcut, lfc){
  nm <- lapply(result, function(r){
    which(r$padj < pcut & abs(r$log2FoldChange) > lfc)
  })
  names(nm) <- lapply(result, getComparison)
  nm
}

#' extract counts of significanlty deregulated genes for one comparison
#' @param comparison deseq2 result tibble
#' @param pcut adjusted p-value cutoff
#' @param lfc log fold change cutoff
#'
#' @export
extractUpDownSignif <- function(comparison, pcut, lfc){
   up <- sum(comparison$padj < pcut & comparison$log2FoldChange > lfc, na.rm = TRUE)
   down <- sum(comparison$padj < pcut & comparison$log2FoldChange < -lfc, na.rm = TRUE)
   dereg <- up + down
   tibble::tibble(p.adj=pcut, log2FC=lfc, dereg=dereg, up=up, down=down)
}

#' extract counts of significanlty deregulated genes for a list of comparisons
#'
#' @export
extractMultiResultsSignif <- function(result, pcut, lfc){
  do.call("rbind", lapply(result, function(r){
      sig <- extractUpDownSignif(r, pcut, lfc)
      cbind(tibble::tibble(comparison=getComparison(r)), sig)
  }))
}

#' extract counts of significanlty deregulated genes for one of comparison
#' but multiple log2 fold changes
#'
#' @export
extractMultiLFCSignif <- function(comparison){
  do.call("rbind", lapply(c(0, 1, 2, 3 ,4), function(lfc){
     extractUpDownSignif(comparison, pcut, lfc)
  }))
}


#' adds labels to chunks in deseq2_pairwise.Rmd (tpl)
#'
#' https://github.com/yihui/knitr-examples/blob/master/041-label-i.Rmd
#'
#'
#' @export
knitLabels <- function(titfile, knitdir){
      knitr::pat_brew() ## switching to pat_brew
      infile <- system.file('deseq2_pairwise.tpl',package="deseq2analysis")
      titfiles <- gsub("_","", titfile)
      withLabels <- knitr::knit(text = readLines(infile))
      outpath <- paste(knitdir,"/comp_",titfile,".Rmd",sep="")
      fileConn<-file(outpath)
      writeLines(withLabels, fileConn)
      close(fileConn)
      knitr::pat_md()
      outpath
}

#' generates a link for md
#'
#' @export
generateMDLink <- function(title, link){
   paste("[",title,"]", "(", link,")",sep="")
}

#' returns a table with output files based on grouping and config
#'
#' @export
outputFilesTable <- function(deseqconfig, deseq.r){
    purrr::map_df(deseq.r, function(deseq.result){
        #comparison was added by me to metadata(deseq.result)
        comparisongroups <- S4Vectors::metadata(deseq.result)$comparison[2:3]
        tit <- getComparison(deseq.result)
        titfile <- getComparisonString(deseq.result)
        resv <- c(comparison=tit)
        if(deseqconfig$output$savetables){
            exl <- paste(titfile, ".diff.norm.xlsx",sep="")
            tab <- paste(titfile, ".diff.norm.tab",sep="")
            exll <- generateMDLink("excel file", exl)
            tabl <- generateMDLink("tab delimited", tab)
            comb <- paste(exll, tabl)
            resv <- c(resv, `DGE tables`=comb)
        }
        if(!is.null(deseqconfig$functional)){
            func <- paste(titfile,"_wt_functional_analysis.html",sep="")
            funcz <-  paste(titfile,"_wt_functional_analysis.zip",sep="")
            funcl <- generateMDLink("html", func)
            funczl <- generateMDLink("zipped tab delimited", funcz)
            comb <- paste(funcl, funczl)
            resv <- c(resv, `functional analysis`=comb)
        }
        dplyr::bind_rows(resv)
    })
}


#' zips all output files
#'
#'
#' @export
zipall <- function(deseqconfig, knitdir){
    zipbase <- gsub(".html",".zip", deseqconfig$outputname)
    zipfile <- paste(knitdir, "/", zipbase ,sep="")
    zips <- dir(knitdir, pattern="*functional_analysis_files.zip", full.names = TRUE)
    html <- dir(knitdir, pattern="*html", full.names = TRUE)
    dn <- dir(knitdir, pattern="*diff.norm*", full.names=TRUE)
    allz <- c(zips, html, dn)
    zip::zipr(zipfile = zipfile, files = allz)
}




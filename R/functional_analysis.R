

#' read species table from inst
#'
#' @export
speciesTable <- function(){
  col_types <- readr::cols(
     species = readr::col_character(),
     database = readr::col_character(),
     common = readr::col_character(),
     clusterprofiler = readr::col_character(),
     msigdbr = readr::col_character()
   )
   readr::read_csv(system.file("species.csv", package="deseq2analysis"), col_types=col_types, progress = FALSE) %>% dplyr::mutate(abase=paste("org.", database, sep=""))
}

#' reads species table and filters based on species name
#'
#' @param speciesName = Mm Hs At Dm ...
#' @param mapdb the map name to extract
#'
#' @export
species2map <- function(speciesName, mapdb){
  m <- speciesTable() %>% dplyr::filter(species == speciesName) %>% dplyr::pull(!!rlang::sym(mapdb))
  if(length(m) == 1){
    m
  } else {
    NA
  }
}


#' map of species to various database and species identifiers
#'
#' @param speciesName short name of species Mm Hs At Dm ..
#'
#' default: "org.speciesName.eg.db"
#'
#' @export
speciesIDTypeToLib <- function(speciesName){
  abase <- species2map(speciesName, "abase")
  clusterprof <- species2map(speciesName, "clusterprofiler")
  msigdbr <- species2map(speciesName, "msigdbr")
  if(is.na(abase)){
    lib <- paste("org.", speciesName, ".eg.db",sep="")
  }else{
    lib <- paste(abase,".db",sep="")
  }
  list(lib=lib,clusterprofiler=clusterprof,msigdbr=msigdbr)
}


#' adds entrez id column to table based on geneid
#'
#' @param resultTable with geneid column
#' @param fromType one of the validIDtypes
#' @param orgdb the database e.g. org.MM.eg.db
#'
#' can return rows with duplicated geneids if multiple entrez for this gene id
#' can return NAs in entrez
#' can return rows with duplicated entrez if multiple geneids for entrez
#'
#' @export
getEntrezIds <- function(resultTable, fromType, libstr){
   library(libstr, character.only = TRUE)
   lib <- get(libstr)
   entrez <- AnnotationDbi::select(lib, keys=resultTable$geneid, columns=c("ENTREZID","GENENAME","SYMBOL"), keytype=fromType)
   tb <- entrez %>% dplyr::group_by_(fromType) %>%
          dplyr::summarize(entrez=ENTREZID[1],entrezConcat=paste(unique(ENTREZID),sep=",",collapse=""),entrezCount=length(unique(ENTREZID)))
                       #   symbol=SYMBOL[1],symbolConcat=paste(unique(SYMBOL),sep=","),symbolCount=length(unique(SYMBOL)),
                      #     genename=GENENAME[1],genenameConcat=paste(unique(GENENAME),sep=","),symbolCount=length(unique(GENENAME)))
   resultTable %>% dplyr::inner_join(tb, by=c(geneid=fromType))
}


#' filters gene universe by no NA entrez, no NA GO mapping
#' meanFilter of deseq2
#'
#' @param table the entrez annotatet table
#' @param libstr the annotation db eg org.Hs.eg.db
#'
#' @export
filterGenesGO <- function(table, libstr){
  library(libstr, character.only = TRUE)
  lib <- get(libstr)
  ids <- table %>% dplyr::filter(!is.na(entrez)) %>% dplyr::pull(entrez)
  wg <- AnnotationDbi::select(lib, keys=ids, columns="GO", keytype="ENTREZID")$ENTREZID
  withGo <- table[table$entrez %in% as.character(wg),]
  withMean <- withGo[!withGo$meanFilter,]
  withMean
}

#' get collapsed GSEA pathways
#'
#' @param gseaL the complete gsea result
#'
#' @export
getGSEACollapsedPathways <- function(gseaL, pvalGo){
  #LOG("called collapse")
  gsea <- gseaL$gsea
  filtered <- gsea %>% dplyr::filter(padj < pvalGo) %>% dplyr::arrange(pval)
  #LOG("filtered collapse")
  collapsed <- fgsea::collapsePathways(filtered, gseaL$pathways, gseaL$stats) #takes long
  #LOG("collapsed")
  collapsed$mainPathways
}

######## GOstats #######

#' tests for GO overrepresentation with GOStats
#'
#'
#' @export
testGO <- function(signifTable, entrezUniverse, ontology, libMap, pvalGo){
  signifEntrez <- signifTable$entrez
  if(length(signifEntrez) > 0){
    tryCatch(
      {
        params <- new("GOHyperGParams", geneIds=signifEntrez, universeGeneIds=entrezUniverse, annotation=libMap$lib,
                      ontology=ontology, pvalueCutoff=pvalGo, conditional=TRUE, testDirection="over")

        hgOver <- hyperGTest(params)
        extractLong <- extract(hgOver, pvalGo, signifTable)
        if(!is.na(extractLong)){
          extractLong %>% dplyr::mutate(ontology=ontology)
        } else {
          NA
        }
      }, error = function(err){ LOG(paste("ERROR but continuing",err)); NA }
    )
  } else {
      NA
  }
}

#' tests for overrrepresentation in KEGG with GOstats
#'
#'
#' @export
testKEGG <- function(signifTable, entrezUniverse, libMap, pvalGo){
  signifEntrez <- signifTable$entrez
  params <- new("KEGGHyperGParams", geneIds=signifEntrez, universeGeneIds=entrezUniverse, annotation=libMap$lib,
                pvalueCutoff=pvalGo, testDirection="over")

  hgOver <- hyperGTest(params)
  kegg <- summary(hgOver) # %>% dplyr::mutate(ontology=ontology)
  kegg
}

#' tests for overrepresentation with ReactomePA
#'
#'
#' @export
testReactome <- function(signifTable, organism, pvalGo){
  signifEntrez <- signifTable$entrez
  ep <- ReactomePA::enrichPathway(gene=signifEntrez, organism=organism, pvalueCutoff=pvalGo, readable=T)
  ep
}


#' tests for functional enrichment GO, kegg, react
#' saves result in RDS
#'
#' @export
testFunctionalEnrich <- function(signifTable, entrezUniverse, libMap, organism, pvalGo, regtype, outname){
  LOG(paste("organism", organism, pvalGo, outname))

  mf <- testGO(signifTable, entrezUniverse, "MF", libMap, pvalGo)
  bp <- testGO(signifTable, entrezUniverse, "BP", libMap, pvalGo)
  cc <- testGO(signifTable, entrezUniverse, "CC", libMap, pvalGo)

  kegg <- NULL
  react <- NULL
  #TODO: detect PATH in keytypes and react ..
  if(FALSE){ # keytypes(libMap$lib) contains PATH ...){
      kegg <- testKEGG(signifTable, entrezUniverse, libMap, pvalGo)
  }
  if(FALSE){#! is.null(organism)){
    react <- testReactome(signifTable, organism, pvalGo)
  }
  li <- list(regtype=regtype,mf=mf,bp=bp,cc=cc,kegg=kegg,react=react)
  saveRDS(li, paste(outname,".RDS",sep=""))
  invisible(li)
}


#' GSEA test
#'
#' saves result in RDS
#'
#' @export
testGSEA <- function(resultTable, organism, pvalGo, outname){
  sortedTable <- resultTable[order(resultTable$log2FoldChange,decreasing = TRUE),]
  geneList <- sortedTable$log2FoldChange
  #LOG("got gene list")
  names(geneList) <- sortedTable$entrez
  geneList <- Filter(function(x){!is.na(x)}, geneList)
  #LOG(paste("filtered gene list",organism))
  geneSetsTable <- msigdbr::msigdbr(species = organism)
  #LOG("extracted geneSets")
  geneSets <- geneSetsTable %>% split(x = .$entrez_gene, f = .$gs_name)
  geneSetsCat <- geneSetsTable %>%
    dplyr::select(gs_name, gs_id, gs_cat) %>%
    dplyr::group_by(gs_name, gs_id, gs_cat) %>%
    dplyr::summarise(cat.size=dplyr::n()) %>%
    dplyr::ungroup()   ## gs_subcat
  #LOG("got catalogs")
  gsea <- fgsea::fgsea(pathways=geneSets, stats=geneList, minSize=15, maxSize=500, nperm=10000)
  gseaj <- gsea %>% dplyr::inner_join(geneSetsCat,by=c(pathway="gs_name"))
  #LOG("joined catalogs")
  gseaL <- list(pathways=geneSets,stats=geneList,gsea=gseaj,pathway_table=geneSetsCat)
  #LOG("before collapse")
  collapsed <- getGSEACollapsedPathways(gseaL, pvalGo) ## takes long
  gseaL$collapsed_pathways <- collapsed
  gseaL$pvalGo <- pvalGo
  saveRDS(gseaL, outname)
}

#' test functional enrichment with config parameters
#'
#' @param resultTable table of ID, etc..
#' @param outbasefun output base of functional annotation, should be valid path
#' @param deseqconfig deseq config
#'
#' @export
testfunctionalFromConfig <- function(resultTable, outbasefun, deseqconfig, comparisontitle){
    functional <- deseqconfig$functional
    organism <- functional$organism
    idtype <- functional$idtype
    pvalExp <- functional$adjp
    pvalGo <- functional$pval
    log2FCExp <- functional$log2fc
    deseqconfig$comparisontitle <- comparisontitle #used for report
    saveRDS(deseqconfig, paste(outbasefun, "_deseq2parameters.RDS", sep=""))
    testFunctional(resultTable, organism, idtype, outbasefun, pvalExp, log2FCExp, pvalGo)
}

#' test functional enrichment
#'
#' @param resultTable table of ID, etc..
#' @param organism the organism abbreviated name: Mm, Hs, Dm, Ce, At
#' @param idtype type of id column in resultTable ENSEMBL etc..
#' @param outbasefun output base of functional annotation, should be valid path
#' @param pvalExp p-value cutoff for expression values
#' @param log2FCExp log2FC for expression values
#' @param pvalGo p-value cutoff for GO analysis
#'
#'
#' @return nothing
#'
#' saves multiple RDS files with GO etct..
#' @import GOstats
#' @import GO.db
#'
#' @export
testFunctional <- function(resultTable, organism, idtype, outbasefun, pvalExp, log2FCExp, pvalGo){
  libMap <- speciesIDTypeToLib(organism)
  libMapStr <- paste(libMap, collapse=" ", paste="")
  #LOG(paste("species annotation map: ", libMapStr))

  #LOG("testFunctional")
  resultEntrez <- getEntrezIds(resultTable, idtype, libMap$lib)
  filteredTable <- filterGenesGO( resultEntrez, libMap$lib )
  entrezUniverse <- unique(filteredTable$entrez)


  signifup <- subset(filteredTable, padj < pvalExp & log2FoldChange >  log2FCExp )
  signifdn <- subset(filteredTable, padj < pvalExp & log2FoldChange < -log2FCExp )
  signifde <- subset(filteredTable, padj < pvalExp & abs(log2FoldChange) > log2FCExp )

  testFunctionalEnrich(signifup, entrezUniverse, libMap, libMap$clusterprofiler, pvalGo, "upregulated", paste(outbasefun,"_up",sep=""))
  testFunctionalEnrich(signifdn, entrezUniverse, libMap, libMap$clusterprofiler, pvalGo, "downregulated", paste(outbasefun,"_dn",sep=""))
  testFunctionalEnrich(signifde, entrezUniverse, libMap, libMap$clusterprofiler, pvalGo, "deregulated", paste(outbasefun,"_de",sep=""))

  #LOG("after signif filter"
  if(FALSE){ #! is.null(libMap$msigdbr)){
    testGSEA(resultEntrez, libMap$msigdbr, pvalGo, paste(outbasefun,"_gse.RDS",sep=""))
  }

}

#' helper function for GO table
#' extracts the results in long format (i.e. one row per geneid that is in category)
#'
#'
#' @export
extract <- function(hyper, maxPval, signifTable, categorySize=NULL){
      wanted <- getWantedResults(hyper, maxPval, categorySize)
      pvals.all <- GOstats::pvalues(hyper)
      ucounts.all <- GOstats::universeCounts(hyper)
      if (!any(wanted)) {
        NA
      } else {
        pvals <- pvals.all[wanted]
        ucounts <- ucounts.all[wanted]
        catIds <- names(pvals)
        odds <- GOstats::oddsRatios(hyper)[wanted]
        ecounts <- GOstats::expectedCounts(hyper)[wanted]
        counts <- GOstats::geneCounts(hyper)[wanted]
        catrows <- lapply(1:length(catIds), function(index){
          goid <- catIds[index]
          #goterm <- AnnotationDbi::Term(goid) ReactomePA hides this
          goterm <- AnnotationDbi:::.GOid2go_termField(goid,"term")
          nd <- graph::nodeData(GOstats::goDag(hyper))[[goid]]
          geneIds <- nd$condGeneIds
          fil <- subset(signifTable, entrez %in% geneIds)
          entrez <- fil$entrez
          gene.ids <- as.character(fil$geneid)
          link <- paste0('<a target="_blank" href="http://amigo.geneontology.org/amigo/term/',goid,'">',goterm,'</a>',sep="")
          tibble::tibble(go.id=goid, go.term=goterm, go.link=link, go.pval=nd$pvalue, go.oddsRatio=nd$oddsRatio,
                             go.expCount=nd$expCount, go.count=counts[index], go.size=ucounts[index],
                             go.gene.id=gene.ids,go.gene.ids=paste(gene.ids,collapse=","),go.result.index=index)
        })
        resultLong <- do.call("rbind", catrows)
        resultLong
      }
}

#' helper function for GO table
#'
#'
#' @export
getWantedResults <- function(result, pvalue, categorySize=NULL){
    ## Returns a logical vector with TRUE indicating selected
    ## results from those tested in the specified result instance.
    pvals <- pvalues(result)
    wanted <- is.finite(pvals) & pvals < pvalue
    if (!is.null(categorySize)) {
      ucounts <- universeCounts(result)
      hasMinSize <- ucounts >= categorySize
      wanted <- wanted & hasMinSize
    }
    wanted
}

#' plots top
#'
#'
#' @export
plotGSEAOverviewTop <- function(gseaL, maxPathway, pvalGo, minSize=20){
  gsea <- gseaL$gsea
  topPathwaysUp <- gsea %>% dplyr::filter(ES > 0 & size > minSize & padj < pvalGo) %>% dplyr::top_n(maxPathway, pval) %>% dplyr::pull(pathway)
  topPathwaysDown <- gsea %>% dplyr::filter(ES < 0 & size > minSize & padj < pvalGo) %>% dplyr::top_n(maxPathway, pval) %>% dplyr::pull(pathway)
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  fgsea::plotGseaTable(gseaL$pathways[topPathways], gseaL$stats, gsea, gseaParam = 0.5)
}


#' extracts collapsed top
#'
#'
#' @export
extractGSEAOverviewCollapsedTop <- function(gseaL, maxPathway, categoryFilter, pvalGo, minSize=20){
  filteredPathways <- gseaL$pathway_table %>% dplyr::filter(categoryFilter == gs_cat & gs_name %in% gseaL$collapsed_pathways) %>% dplyr::pull(gs_name)
  gsea <- gseaL$gsea
  filtered <- gsea %>% dplyr::filter(padj < pvalGo) %>% dplyr::arrange(pval)
  topPathwaysUp <- gsea %>% dplyr::filter(ES > 0 & size > minSize & padj < pvalGo & pathway %in% filteredPathways) %>%
    dplyr::top_n(maxPathway, pval) %>% dplyr::pull(pathway)
  topPathwaysDown <- gsea %>% dplyr::filter(ES < 0 & size > minSize & padj < pvalGo & pathway %in% filteredPathways) %>%
    dplyr::top_n(maxPathway, pval) %>% dplyr::pull(pathway)
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  topPathways
}


#' path of output table
#'
#'
#' @export
outTableFunctionalPath <- function(outtabledir, outbasefun, suffix){
    paste(outtabledir, "/", outbasefun, "_", suffix, "_functional_analysis.tab",sep="")
}


#' scientific format for HTML tables
#'
#' @export
sfformat <- function(value){
   formatC(value, "E", width=3, digits=2)
}

#' datatable for GO
#' data comes in in long format, which is why distinct is necessary
#'
#' also creates a text file
#'
#' @export
datatabGO <- function(dat, suffix, outtabledir, outbasefun){
  tabvalshtml <- dat %>% dplyr::ungroup() %>% dplyr::select(-dplyr::one_of(c("go.result.index","go.oddsRatio","go.term","go.expCount","go.gene.id","go.gene.ids", "ontology"))) %>% dplyr::distinct()
  if(! is.null(outtabledir) & ! is.null(outbasefun)){
     tabvalstab <- dat %>% dplyr::ungroup() %>% dplyr::select(-dplyr::one_of(c("go.result.index", "go.link","go.gene.id")))  %>% dplyr::distinct()
     readr::write_tsv(tabvalstab, outTableFunctionalPath(outtabledir, outbasefun, suffix))
  }
  tabvalshtml %>% dplyr::mutate(pvalue=sfformat(go.pval)) %>% dplyr::select(-go.pval) %>% DT::datatable(escape=FALSE)
}

#' datatable for KEGG
#'
#' also creates a text file
#'
#' @export
datatabKEGG <- function(dat, suffix, outtabledir, outbasefun){
  if(! is.null(outtabledir) & ! is.null(outbasefun)){
     readr::write_tsv(dat, outTableFunctionalPath(outtabledir, outbasefun, suffix))
  }
  dat %>% dplyr::mutate(pvalue=sfformat(Pvalue), `expected count`=sfformat(ExpCount))  %>% dplyr::select(-c(Pvalue, ExpCount)) %>% DT::datatable(escape=FALSE)

}

#' datatable for React
#'
#' also creates a text file
#'
#' @export
datatabReact <- function(dat, suffix, outtabledir, outbasefun){
  if(! is.null(outtabledir) & ! is.null(outbasefun)){
    readr::write_tsv(dat, outTableFunctionalPath(outtabledir, outbasefun, suffix))
  }
  dat %>% dplyr::mutate(p.adj=sfformat(p.adjust), qvalue=sfformat(qvalue), pvalue=sfformat(pvalue)) %>% dplyr::select(-p.adjust) %>% DT::datatable(escape=FALSE)
}

#' datatable for GSEA
#'
#' also creates a text file
#'
#' @export
datatabGSEA <- function(gsea, toppathways, suffix, outtabledir, outbasefun){
  dat <- gsea %>% dplyr::filter(pathway %in% toppathways) %>% dplyr::select("pathway","padj","size","gs_id") #ES,NES
  dats <- dat[match(toppathways, dat$pathway),] %>%
    dplyr::mutate(pathwayl=paste0('<a target="_blank" href="http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=',pathway,'">',pathway,'</a>',sep="")) %>%
    dplyr::select(-pathway)
  if(! is.null(outtabledir) & ! is.null(outbasefun)){
    readr::write_tsv(dats, outTableFunctionalPath(outtabledir, outbasefun, suffix))
  }
  dat  %>% dplyr::mutate(p.adj=sfformat(padj))  %>% dplyr::select(-padj) %>% DT::datatable(escape=FALSE)
}


#' calls function if item is not NA or NULL
#' otherwise returns NO DATA table
#'
#'
#' @export
emptyTable <- function(item, fun, ...){
  if(is.object(item)){
    fun(item, ...)
  }else{
    DT::datatable(data.frame(no=NA,data=NA), escape=FALSE)
  }
}

#' plots react only if data frome has rows
#'
#' also creates a text file
#'
#' @export
plotEmpty <- function(react, pl){
  empty <- nrow(as.data.frame(react)) == 0
  if(empty){
    ggplot(data.frame(text="no data"), aes(x=1,y=1,label=text)) + geom_text() + theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.text=element_blank())
  } else {
    pl
  }
}

#' describes the parameters for the functional analysis
#'
#' @export
goparameterdescription <- function(deseqconfig){
  if(is.null(deseqconfig)){
    ""
  }else{
    functional <- deseqconfig$functional
    organism <- functional$organism
    idtype <- functional$idtype
    pvalExp <- functional$adjp
    pvalGo <- functional$pval
    log2FCExp <- functional$log2fc


    paste("For determining up, down and deregulated genes we used an adjusted p-value cutoff of ", pvalExp, " and a minimal log2 fold change of ", log2FCExp,". ",
          "The p-value cutoff for the GO-Term and KEGG analysis was ",pvalGo,". ", sep="")
  }
}


#' zips all files in the functional_analysis folder
#'
#'
#' @export
zipfunctabs <- function(outtabledir, knitdir){
  if(!is.null(outtabledir)){
        zipfile <- paste(knitdir, "/", outtabledir,".zip",sep="")
        files2zip <- dir(paste(knitdir, "/", outtabledir, sep=""), full.names = TRUE)
        zip::zipr(zipfile = zipfile, files = files2zip)
  }
}


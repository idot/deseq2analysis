


########## plots for alignment stats ############

#' Plots the alignment statistics
#'
#' @export
plotAlignments <- function(alignments){
  alignFiltered <- subset(alignments, countSum > 0)
  p1 <-  ggplot(alignFiltered, aes(x=sampleId, y=countSum, fill=V1)) + geom_bar(stat="identity", position="stack") + ylab("absolute counts") + xlab("sample id") + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + guides(fill=guide_legend(title="align type",ncol=2)) + idoplots::discrete_fill()
  p2 <-  ggplot(alignFiltered, aes(x=V1, y=percent, colour=sampleId)) + geom_point() + ylab("percent of total") + xlab("align type") + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + idoplots::discrete_colour()
  if(length(unique( alignFiltered$sampleId ) > 10)){
    p2 <- p2 + guides(colour=guide_legend(ncol=2))
  }
  gridExtra::grid.arrange(p1, p2, ncol=1)
}

#' Plots the zero
#'
#' @export
plotZero <- function(xplicates){
  p1 <- ggplot(xplicates, aes(x=freq, y=cumPerc, colour=sampleId)) + geom_line() + scale_x_log10(limits=c(1,1e3), breaks=c(1,2,3,5,10,20,100,1000)) + xlab("f
                                                                                                                                                           requency of read position") + ylab("cumulative percentage of total") + idoplots::discrete_colour()
  p2 <- ggplot(subset(xplicates, freq <= 10), aes(x=freq, y=freqCountNorm,colour=sampleId)) + geom_line() + xlab("frequency of read position") + ylab("percentage of total") + scale_x_continuous(breaks=1:10) + idoplots::discrete_colour()
  idoplots::grid_arrange_shared_legend(p1, p2, nrow=2, ncol=1)
}

#' Plots the coverages
#'
#' @export
plotCoverages <- function(coverageLong){
  maxLong <- ceiling(max(coverageLong$mean))
  ggplot(coverageLong, aes(x=bin,y=mean,colour=sampleId)) + geom_line() + facet_grid(. ~ sense) + xlab("position in length normalised mRNA") + ylab("normalized mean coverage") + coord_cartesian(ylim=c(0,maxLong)) + theme(legend.position="bottom") + idoplots::discrete_colour()
}


#' Plots the annotation distribution
#'
#' @export
plotAnnotationDistribution <- function(annot){
  ggplot(annot, aes(x=sampleId, y=countFraction, fill=Group)) + geom_bar(position="stack",stat="identity") + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + ylab("fraction") + idoplots::discrete_fill()
}

#' Plots the dupfit fit
#'
#' @export
plotDupFit <- function(dupfit){
  ggplot(subset(dupfit,fit=="RPKM"), aes(x=intercept,y=slope,colour=sampleId)) + geom_point() + xlab("percentage of duplication at low RPKM (fitted)") + ylab("Slope (log of duprate growth for 1 RPKM)") + idoplots::discrete_colour()
}

########## plots for aggreagate data ####################

#### DEVis TODO
####plot_mds(filename="MDS.pdf", color_var="treatment", shape_var="time", theme=1)




########## plots for differential expression ############


#' vulcano plot
#'
#' converts infinite log2FC values to max(log2FC) * 0.1
#'
#' @export
plotVulcano <- function(comp, pcut, lfc, compname){
  comps <- comp %>% dplyr::filter(is.finite(log2FoldChange)) %>%
          dplyr::mutate( significant = is.finite(padj) & padj < pcut & abs(log2FoldChange) > lfc)
  mvallf <- max(abs(comps$log2FoldChange))
  mvallf <- mvallf * 1.01
  comps_top <- comps %>% dplyr::filter(significant) %>% dplyr::arrange(padj) %>% dplyr::slice(1:10)
  p <- ggplot(comps, aes(x=log2FoldChange,y=-log10(padj),colour=significant)) + geom_point(aes(alpha=ifelse(significant, 0.7,0.3))) +
    geom_text(data=comps_top, aes(x=log2FoldChange,y=-log10(padj),label=geneid),position = position_nudge(x = 0.3,y=0.5), show.legend=FALSE) +
    xlab(compname) + ylab("-log10(p-adjusted)") + scale_alpha_identity() + guides(colour=guide_legend("significant"), alpha=FALSE)
  p
}


#' MA plot
#'
#'
#'
#' @export
plotMAPlot <- function(comp, pcut, lfc, compname){
  maxLFC <- max(abs(comp$log2FoldChange),na.rm=TRUE) * 1.05
  maxLFCs <- maxLFC * 1.01
  comps <- comp %>%
    dplyr::mutate(lfci=!is.finite(log2FoldChange), lfcm = ifelse(lfci, maxLFC * sign(log2FoldChange),log2FoldChange) , significant = is.finite(padj) & padj < pcut & abs(log2FoldChange) > lfc)
  ggplot(comps, aes(x=baseMean,y=lfcm,colour=significant)) +
     geom_point(aes(size=2,alpha=ifelse(significant, 0.7,0.3))) +
     idoplots::scale_x_log2() +
     scale_y_continuous(limits=c(-maxLFCs,maxLFCs)) + xlab("mean of normalized counts") + ylab(compname) +
     scale_alpha_identity() + scale_size_identity() + guides(colour=guide_legend("significant"), alpha=FALSE, size=FALSE)
}

#' p-value distribution plot
#'
#'
#'
#' @export
plotPvalDist <- function(comp, pcut){
  pvaldf <- data.frame(pvalue=c(rep("raw",nrow(comp)), rep("adjusted",nrow(comp))),value=c(comp$pvalue,comp$padj))
  ggplot(pvaldf, aes(x=value, fill=pvalue)) + geom_histogram(binwidth=0.01, position="identity") + facet_grid(pvalue ~ .) + geom_vline(xintercept = pcut, color="red")
}

#' independent filtering plot
#'
#'
#'
#' @export
plotIndependentFiltering <- function(comp, pcut, lfc, treshold){
  ggplot(comp %>% dplyr::mutate(significant=is.finite(padj) & padj < pcut & abs(log2FoldChange) > lfc), aes(x=baseMean, y=-log10(padj),color=significant)) +
           geom_point() +
           geom_vline(xintercept = treshold, color="red") + idoplots::scale_x_log2() +
           ylab("-log10(adjusted p-value)") + xlab("log2(base mean)")
}

#' size factors plot
#'
#'
#'
#' @export
plotSizeFactors <- function(dds.r){
   ggplot(as.data.frame(SummarizedExperiment::colData(dds.r)), aes(x=sampleId, y=sizeFactor,color=group)) + geom_point() + geom_hline(yintercept = 1, color="darkgray") + ylab("size factor") + idoplots::xrot()
}

#' total counts plot
#'
#'
#'
#' @export
plotTotalCounts <- function(dds){
    totalCounts <- colSums(BiocGenerics::counts(dds))
    counts <- tibble::enframe(totalCounts, "sampleId", "count") %>%
         dplyr::left_join(as.data.frame(SummarizedExperiment::colData(dds)), by=c(sampleId="sampleId"))
    ggplot(counts,aes(x=sampleId,y=count,color=group)) + geom_point() + idoplots::xrot()
}


#' plot counts per gene normalised
#'
#'
#' @export
plotCountsPerGene <- function(dds.r, plus=0){
    norm <-  BiocGenerics::counts(dds.r, normalized=TRUE)
    normm <- reshape2::melt(norm, value.name="normalised", varnames = c("geneid", "sampleId")) %>%
        dplyr::mutate(sampleId=as.character(sampleId)) %>%
        dplyr::left_join(as.data.frame(SummarizedExperiment::colData(dds.r)), by=c(sampleId="sampleId"))

    p <- ggplot(normm, aes(x=normalised + plus, color=group, linetype=as.factor(internal_replicate))) + geom_density()
    p  + idoplots::scale_x_log2()
}

#' plot counts per gene normalised (+/- 1)
#'
#'
#' @export
plotCountsPerGeneBoth <- function(dds.r){
    plus <- plotCountsPerGene(dds.r, 1) + xlab("normalised + 1")
    zero <- plotCountsPerGene(dds.r, 0) + xlab("normalised, no 0 counts")
    idoplots::grid_arrange_shared_legend(plus, zero)
}

#' gets data frame with dispersion estimates
#'
#'
#' @export
dispEstToDF <- function(dds.r){
  px <- SummarizedExperiment::mcols(dds.r)$baseMean
  sel <- (px > 0)
  means <- px[sel]

  #to make it comparable across reports same axis
  dispEstXMin <- 1e-1
  dispEstXMax <- 1e6
  dispEstYMin <- 1e-5
  dispEstYMax <- 1e2

  toMinMax <- function(vs, vmin, vmax){
    vs[vs > vmax] <- vmax*0.99
    vs[vs < vmin] <- vmin*2
    vs
  }

  py <- SummarizedExperiment::mcols(dds.r)$dispGeneEst[sel]
  outlier <-SummarizedExperiment::mcols(dds.r)$dispOutlier[sel]
  final <- DESeq2::dispersions(dds.r)[sel]
  fitted <- SummarizedExperiment::mcols(dds.r)$dispFit[sel]
  below <- 2*dispEstYMin ## arbitrary value for very low dispersion genes
  belowT <- py < below
  py[belowT] <- below
  pointsDF <- data.frame(mean=toMinMax(means, dispEstXMin, dispEstXMax),gene=toMinMax(py,dispEstYMin,dispEstYMax),outlier=outlier,fitted=toMinMax(fitted,dispEstYMin,dispEstYMax),final=toMinMax(final,dispEstYMin,dispEstYMax))
  pointsDF
}

#' plot the dispersion estimate
#'
#'
#' @export
plotDispEst <- function(dds){
  #to make it comparable across reports same axis
  dispEstXMin <- 1e-1
  dispEstXMax <- 1e6
  dispEstYMin <- 1e-5
  dispEstYMax <- 1e2

  dispDF <- dispEstToDF(dds)
  p <- ggplot(dispDF, aes(x=mean,y=gene)) + geom_point(size=2,alpha=0.3,aes(colour="gene estimate")) + geom_point(size=2,alpha=0.3,aes(x=mean,y=final,colour="final")) + geom_point(size=1,aes(x=mean,y=fitted,colour="fitted")) + scale_x_log10(limits=c(dispEstXMin,dispEstXMax)) + scale_y_log10(limits=c(dispEstYMin,dispEstYMax)) + guides(colour=guide_legend("category")) + xlab("mean of normalized counts") + ylab("dispersion")
  if(any(dispDF$outlier)){
    p <- p + geom_point(data=subset(dispDF, outlier), alpha=0.4, size=3, aes(x=mean,y=gene,colour="outlier"))
  }
  p
}



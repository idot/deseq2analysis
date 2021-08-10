#!/usr/bin/env Rscript

#R packaging does not allow executable scripts
#this would be a runner script

option_list <- list(
  optparse::make_option(c("-c","--config"), type="character", help="yaml config file"),
  optparse::make_option(c("-k","--knitdir"), type="character", default=".", help="knitting and output folder, very important for singularity"),
  optparse::make_option(c("--keep_intermediate"), action="store_true", default=FALSE, help="keep intermediate knit documents for debugging")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list),positional_arguments =FALSE)
deseq2analysis::analyseFromConfig(opt$config, opt$knitdir, opt$keep_intermediate)





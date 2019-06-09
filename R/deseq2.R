
#' runs the whole analysis from a yaml config
#' see files in inst/test
#'
#' @param configpath path to yaml config
#' @param knitdir folder to knit and create output in (caveat singularity. wdir is not writable)
#'
#' @export
analyseFromConfig <- function(configpath, knitdir, keep_intermediate=FALSE){
		deseqconfig <- yaml::read_yaml(configpath)
		KNITDIR <- knitdir
		CLEAN <- !keep_intermediate
    print(deseqconfig)

		form <- bookdown::gitbook( self_contained = TRUE, split_by="none", config = list() )

		rmarkdown::render(system.file("deseq2.Rmd", package="deseq2analysis"), knit_root_dir=KNITDIR, intermediates_dir=KNITDIR,
                 output_format=form, output_file=deseqconfig$outputname, output_dir=KNITDIR, clean=CLEAN,
                 params = list(analysis_title = deseqconfig$analysis_title )  )


		rdsfiles <- gsub("_deseq2parameters.RDS", "", dir(".","*_deseq2parameters.RDS"))
		for(outbasefun in rdsfiles){
		      paramfile <- paste(outbasefun,"_deseq2parameters.RDS",sep="")
		      params <- readRDS(paramfile)
    			outbasehtmlfinal <- paste(outbasefun,"_functional_analysis.html",sep="")
    				rmarkdown::render(system.file("functional_analysis.Rmd", package="deseq2analysis"), knit_root_dir=KNITDIR, intermediates_dir=KNITDIR,
                	output_format=form, output_file=outbasehtmlfinal, output_dir=KNITDIR, clean=CLEAN, envir = new.env(),
                	params = list(functional_title = paste(params$comparisontitle, "functional analysis"))
                	)
		}


}



.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis/lib/R/library"))
.libPaths()
# sessionInfo()
Sys.setenv(RETICULATE_PYTHON="/home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis/bin/python")

suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(dplyr)) 
suppressMessages(library(ggplot2)) 
suppressMessages(library(Matrix)) 
suppressMessages(library(data.table)) 
suppressMessages(library(ggpubr)) 
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library("optparse"))
suppressMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(future))


opt = NULL

options(warn = 1)

log_info_simple <- function(message) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(timestamp, "INFO:", message, "\n")
}


data_wrangling = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  

  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ---- "DE_results_without_time_as_a_covariate_NEW_METHOD_only_padj_less_0.05.tsv"
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform processors ----
  
  processors = opt$processors
  
  cat("processors\n")
  cat(sprintf(as.character(processors)))
  cat("\n")
  
  #### READ and transform memory ----
  
  #### READ and transform total_memory (memory in MB) ----
  total_memory = opt$total_memory # Corrected variable name to match bash script
  cat("Total Memory (MB) for global objects:", as.character(total_memory), "\n") # Improved log message
  
  #### Assign resources -------------
  
  log_info_simple("plan stage")
  
  # Set up parallel processing: 'multiprocess' works on a single machine across cores.
  # 'total_memory' is expected in MB from the bash script, convert to bytes for future.globals.maxSize.
  plan("multicore", workers = processors)
  options(future.globals.maxSize = total_memory * 1024^2) # Corrected: Convert MB to bytes
 
  #### READ adata ----
  
  adata<-readRDS(file=opt$SeuratObject)
  

  # cat("adata_0\n")
  # cat(str(adata))
  # cat("\n")
  
  #### READ Linked_peak_to_DE_genes_ATAC_by_Construction_annotation ----
  
  Linked_peak_to_DE_genes_ATAC_by_Construction_annotation<-readRDS(file=opt$Linked_peak_to_DE_genes_ATAC_by_Construction_annotation)
  
  
  cat("Linked_peak_to_DE_genes_ATAC_by_Construction_annotation\n")
  cat(str(Linked_peak_to_DE_genes_ATAC_by_Construction_annotation))
  cat("\n")
 
  # Set ATAC_by_Construction_annotation as assay ----------------------
 
  DefaultAssay(adata)<-'ATAC_by_Construction_annotation'
  
  peaks <- GetAssayData(adata,slot="counts")
  
  Peak_list=data.frame(rownames(peaks))
  Peak_list$peak_ids=Peak_list$rownames.peaks.
  
  Peak_list_filt=Peak_list %>% dplyr::filter(peak_ids %in% Linked_peak_to_DE_genes_ATAC_by_Construction_annotation$Peak_ID)
  
  cat("Peak_list_filt_0\n")
  cat(str(Peak_list_filt))
  cat("\n")
  
  peak_subset = peaks[Peak_list_filt$peak_ids,]
  
  cat("peak_subset_0\n")
  cat(str(peak_subset))
  cat("\n")
  
  adata[["peak_subset"]] <- CreateChromatinAssay(
    counts = peak_subset,
    fragments = Fragments(adata),
    min.cells = 0, min.features = -1 )
  
  
  ########### save  --------------
  
  setwd(out)
  
  saveRDS(adata,file="merged_clusters_final_with_peak_subset.rds")
  
  

}


printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}



#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--SeuratObject"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Linked_peak_to_DE_genes_ATAC_by_Construction_annotation"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--processors"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--total_memory"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  data_wrangling(opt)
 

}


###########################################################################

system.time( main() )
  

.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/lib/R/library"))
.libPaths()
# sessionInfo()

Sys.setenv(RETICULATE_PYTHON="/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/bin/python")
library(reticulate)
reticulate::use_python("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/bin/python")
reticulate::use_condaenv("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF")
reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg')
suppressMessages(library("optparse"))
suppressMessages(library(hdf5r))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(scDblFinder))
suppressMessages(library("tidyr"))
suppressMessages(library("tibble"))
suppressMessages(library("biovizBase"))
suppressMessages(library("patchwork"))
suppressMessages(library(glmGamPoi))
suppressMessages(library(future))


log_info_simple <- function(message) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(timestamp, "INFO:", message, "\n")
}





opt = NULL

options(warn = -1)

merge_and_recall_peaks = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform processors ----
  
  processors = as.numeric(opt$processors)
  
  cat("processors\n")
  cat(sprintf(as.character(processors)))
  cat("\n")
  
  #### READ and transform memory ----
  
  #### READ and transform total_memory (memory in MB) ----
  total_memory = as.numeric(opt$total_memory) # Corrected variable name to match bash script
  cat("Total Memory (MB) for global objects:", as.character(total_memory), "\n") # Improved log message
  
  #### Assign resources -------------
  
  log_info_simple("plan stage")
  
  # Set up parallel processing: 'multiprocess' works on a single machine across cores.
  # 'total_memory' is expected in MB from the bash script, convert to bytes for future.globals.maxSize.
  plan("multicore", workers = processors)
  options(future.globals.maxSize = total_memory * 1024^2) # Corrected: Convert MB to bytes
  
  # Define cluster column name *before* the first save, using res_param
  cluster_col_name = paste0("RNA_snn_res.", opt$res_param)
 
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Read filtered object by doublets -----
  
  
  adata2<-readRDS(file=opt$db_filt_clustered_QCed)
  
  cat("adata2_0\n")
  cat(str(adata2))
  cat("\n")
  
  
  
  #### Cluster without integration and check if any cluster looks specially bad for QC metrics -------------
  
  DefaultAssay(adata2) <- 'RNA'
  
  adata2 <- SCTransform(adata2, verbose = FALSE) 
  adata2 <- RunPCA(adata2) 
  adata2 <- RunUMAP(adata2, dims=1:50, reduction.name='umap.rna', reduction.key='rnaUMAP_')
  
  #### ATAC modality -------------
  # We exclude the first dimension as this is typically correlated with sequencing depth
  
  
  DefaultAssay(adata2) <- 'ATAC'

  adata2 <- RunTFIDF(adata2)
  adata2 <- FindTopFeatures(adata2, min.cutoff='q0')
  adata2 <- RunSVD(adata2)
  adata2 <- RunUMAP(adata2, reduction='lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')
  
  
  #### WNN ATAC+RNA modality -------------
  
  adata2 <- FindMultiModalNeighbors(adata2, reduction.list=list('pca', 'lsi'), dims.list=list(1:50, 2:50))
  adata2 <- RunUMAP(adata2, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
  adata2 <- FindClusters(adata2, graph.name='wsnn', algorithm=4, resolution = opt$res_param, verbose=FALSE, method = "igraph")
  
  
  ###### SAVE -----
  
  setwd(out)
  
  saveRDS(adata2, file = 'merged_unprocessed_db_filt_clustered_QCed_reclustered.rds')
  
  # --- DIAGNOSTIC AND PREPARATION ---
  
  # 1. Re-extract the matrix (using the successful logic from the last attempt)
  assays_list <- names(adata2)
  assay_name <- ifelse("SCT" %in% assays_list, "SCT", "RNA")
  slot_name <- ifelse(assay_name == "SCT", "data", "scale.data")
  
  corrected_sct_matrix <- GetAssayData(adata2, assay = assay_name, slot = slot_name)
  
  # 2. Force the matrix to the standard dgCMatrix format (if not already)
  if (!inherits(corrected_sct_matrix, "dgCMatrix")) {
    corrected_sct_matrix <- as(corrected_sct_matrix, "dgCMatrix")
  }
  
  # 3. CRITICAL DIAGNOSTIC: Check for Non-Finite values (NaN, Inf)
  #    These will cause the 'unknown type' HDF5 error.
  num_nan <- sum(is.nan(corrected_sct_matrix@x))
  num_inf <- sum(is.infinite(corrected_sct_matrix@x))
  
  if (num_nan > 0 || num_inf > 0) {
    cat(paste0("WARNING: Found ", num_nan, " NaN values and ", num_inf, " Inf values. These MUST be removed.\n"))
    
    # Replace non-finite values with 0 (or a very small number, depending on your downstream needs)
    # Replacing with 0 is standard for sparse matrices.
    corrected_sct_matrix@x[is.nan(corrected_sct_matrix@x)] <- 0
    corrected_sct_matrix@x[is.infinite(corrected_sct_matrix@x)] <- 0
  } else {
    cat("Matrix passed non-finite value check. Proceeding with export.\n")
  }
  
  # --- R CODE: Save the cleaned objects to disk ---
  
  # 1. Finalize the Seurat object to ensure all components are aligned
  assays_list <- names(adata2)
  assay_name <- ifelse("SCT" %in% assays_list, "SCT", "RNA")
  slot_name <- ifelse(assay_name == "SCT", "data", "scale.data")
  
  corrected_sct_matrix <- GetAssayData(adata2, assay = assay_name, slot = slot_name)
  if (!inherits(corrected_sct_matrix, "dgCMatrix")) {
    corrected_sct_matrix <- as(corrected_sct_matrix, "dgCMatrix")
  }
  corrected_sct_matrix@x <- as.numeric(corrected_sct_matrix@x)
  
  cell_barcodes_to_keep <- colnames(corrected_sct_matrix)
  cell_metadata <- adata2@meta.data[cell_barcodes_to_keep, ] 
  
  # 2. Save the matrix and metadata separately
  setwd(opt$out)
  
  saveRDS(corrected_sct_matrix, file = "final_sct_matrix.rds")
  saveRDS(cell_metadata, file = "final_cell_metadata.rds")
  
  cat("SUCCESS: Saved matrix and metadata as separate RDS files. Proceed to Python.\n")
  
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
    make_option(c("--db_filt_clustered_QCed"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--res_param"), type="numeric", default=NULL, 
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
  
  merge_and_recall_peaks(opt)
 

}


###########################################################################

system.time( main() )
  
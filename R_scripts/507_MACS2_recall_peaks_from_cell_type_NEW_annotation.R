
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

MACS2_call_peaks = function(option_list)
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
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
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
  
 
  #### Read db_filt_clustered_QCed_cell_annotated_MJ_only_genotyped_reclustered -----
  
  
  adata<-readRDS(file=opt$db_filt_clustered_QCed_cell_annotated_MJ_only_genotyped_reclustered)
  
  # cat("adata_0\n")
  # cat(str(adata))
  # cat("\n")



  ## Call Peaks and make new peak matrix


  DefaultAssay(adata) <- 'ATAC'
  peaks <- CallPeaks(
    object = adata,
    group.by = "Construction_annotation",
    macs2.path = "/group/soranzo/conda_envs/Manuel_macs2/bin/macs2")

  frag_file<-opt$frag_file

  Fragmobj <- CreateFragmentObject(frag_file,cells =Cells(adata))


  peakmat = FeatureMatrix(fragments = Fragmobj, features = peaks, cells = Cells(adata), process_n = 10000,
                          sep = c(":", "-"), verbose = TRUE)

  norm_chr = rownames(peakmat)[stringr::str_split_fixed(rownames(peakmat), ":",2)[,1] %in%
                                 paste0("chr", c(1:22, "X", "Y"))]

  peakmat=peakmat[norm_chr,]

  suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
  seqlevelsStyle(annotations)  <- 'UCSC'
  genome(annotations)          <- 'hg38'

  suppressWarnings(chrom_assay <- CreateChromatinAssay(counts=peakmat, sep=c(':', '-'),
                                                       genome='hg38', fragments=Fragmobj,
                                                       min.cells=-1, min.features=-1,
                                                       annotation=annotations))



  adata[['ATAC_by_Construction_annotation']] <- chrom_assay

  #### Save third clustered filtered object  ------------------


  setwd(out)

  saveRDS(adata,file="merged_clusters_final.rds")
  
  # --- DIAGNOSTIC AND EXPORT RNA/SCT MATRIX ---
  
  log_info_simple("Processing RNA/SCT matrix for export.")
  assays_list <- names(adata)
  
  # V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V 
  # CRITICAL MODIFICATION: Use the 'RNA' assay and the 'counts' slot to export unnormalized data
  rna_assay_name <- "RNA"
  rna_slot_name <- "counts" # Exporting RAW COUNTS (Cell Bender corrected)
  # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
  
  corrected_rna_matrix <- GetAssayData(adata, assay = rna_assay_name, slot = rna_slot_name)
  
  if (!inherits(corrected_rna_matrix, "dgCMatrix")) {
    corrected_rna_matrix <- as(corrected_rna_matrix, "dgCMatrix")
  }
  
  # QC: Check for Non-Finite values and replace with 0
  num_nan_rna <- sum(is.nan(corrected_rna_matrix@x))
  num_inf_rna <- sum(is.infinite(corrected_rna_matrix@x))
  
  if (num_nan_rna > 0 || num_inf_rna > 0) {
    cat(paste0("WARNING: Found ", num_nan_rna, " NaN values and ", num_inf_rna, " Inf values in RNA matrix. Removing them.\n"))
    corrected_rna_matrix@x[is.nan(corrected_rna_matrix@x)] <- 0
    corrected_rna_matrix@x[is.infinite(corrected_rna_matrix@x)] <- 0
  } else {
    cat("RNA/SCT matrix passed non-finite value check.\n")
  }
  
  corrected_rna_matrix@x <- as.numeric(corrected_rna_matrix@x)
  
  # 1. Save the RNA matrix
  saveRDS(corrected_rna_matrix, file = "final_rna_corrected_unormalized_matrix.rds")
  log_info_simple("Saved final_rna_corrected_unormalized_matrix.rds")
  
  quit(status = 1)
  
  # --- NEW CODE: DIAGNOSTIC AND EXPORT ATAC MATRIX ---
  
  log_info_simple("Processing ATAC_by_Construction_annotation matrix for export.")
  
  # Exporting the processed/normalized ATAC data (TFIDF-normalized)
  atac_matrix <- GetAssayData(adata, assay = "ATAC_by_Construction_annotation", slot = "data") 
  
  if (!inherits(atac_matrix, "dgCMatrix")) {
    atac_matrix <- as(atac_matrix, "dgCMatrix")
  }
  
  # QC: Check for Non-Finite values in ATAC matrix
  num_nan_atac <- sum(is.nan(atac_matrix@x))
  num_inf_atac <- sum(is.infinite(atac_matrix@x))
  
  if (num_nan_atac > 0 || num_inf_atac > 0) {
    cat(paste0("WARNING: Found ", num_nan_atac, " NaN values and ", num_inf_atac, " Inf values in ATAC matrix. Removing them.\n"))
    atac_matrix@x[is.nan(atac_matrix@x)] <- 0
    atac_matrix@x[is.infinite(atac_matrix@x)] <- 0
  } else {
    cat("ATAC matrix passed non-finite value check.\n")
  }
  
  atac_matrix@x <- as.numeric(atac_matrix@x)
  
  # --- FINAL EXPORT ---
  
  cell_barcodes_to_keep <- colnames(corrected_rna_matrix) # Use RNA barcodes for consistent order
  cell_metadata <- adata@meta.data[cell_barcodes_to_keep, ] 
  
 
  
  # 2. Save the ATAC matrix (new file)
  saveRDS(atac_matrix, file = "final_atac_matrix.rds")
  log_info_simple("Saved final_atac_matrix.rds.")
  
  # 3. Save the metadata
  saveRDS(cell_metadata, file = "final_cell_metadata.rds")
  log_info_simple("Saved final_cell_metadata.rds.")
  
  cat("SUCCESS: Saved RNA/SCT, ATAC matrices, and metadata as separate RDS files. Proceed to Python.\n")
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
    make_option(c("--db_filt_clustered_QCed_cell_annotated_MJ_only_genotyped_reclustered"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--frag_file"), type="character", default=NULL, 
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
  
  MACS2_call_peaks(opt)

 

}


###########################################################################

system.time( main() )
  
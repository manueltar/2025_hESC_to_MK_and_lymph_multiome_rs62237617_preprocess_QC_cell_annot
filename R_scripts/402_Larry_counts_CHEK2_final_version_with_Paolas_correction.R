
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggupset", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("glmnet", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("pheatmap", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("stringr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggsci",lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("patchwork", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))



opt = NULL

options(warn = 1)

read_all_barcoded_samples = function(option_list)
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
  
  #### READ samples ----
  
  samples = unlist(strsplit(opt$samples, split=","))
  
  cat("samples_\n")
  cat(sprintf(as.character(samples)))
  cat("\n")
  
  #### READ genotypes_string ----
  
  genotypes_string = unlist(strsplit(opt$genotypes_string, split=","))
  
  cat("genotypes_string_\n")
  cat(sprintf(as.character(genotypes_string)))
  cat("\n")
  
  #### READ and transform indir ----
  
  indir = opt$indir
  
  cat("indir_\n")
  cat(sprintf(as.character(indir)))
  cat("\n")
  
  #### READ and transform Threshold_attributed_genotypes ----
  
  Threshold_attributed_genotypes = opt$Threshold_attributed_genotypes
  
  cat("Threshold_attributed_genotypes_\n")
  cat(sprintf(as.character(Threshold_attributed_genotypes)))
  cat("\n")
  
  #### READ and transform Threshold_UMIS_per_cell ----
  
  Threshold_UMIS_per_cell = opt$Threshold_UMIS_per_cell
  
  cat("Threshold_UMIS_per_cell_\n")
  cat(sprintf(as.character(Threshold_UMIS_per_cell)))
  cat("\n")
 
  
  #### Create the data frame to read all the file ----
  
  setwd(indir)
  
  files_df<-list.files()
  
  cat("files_df_0\n")
  cat(str(files_df))
  cat("\n")
  
  files_df_sel<-files_df[grep("\\.all_geno_bc_umi$",files_df)]
  
  cat("files_df_sel_0\n")
  cat(str(files_df_sel))
  cat("\n")
  
  samples_files_df_sel<-gsub("\\.all_geno_bc_umi$","",files_df_sel)
  
  cat("samples_files_df_sel_0\n")
  cat(str(samples_files_df_sel))
  cat("\n")
  
  
  
  df_files <- as.data.frame(cbind(samples_files_df_sel,files_df_sel), stringsAsFactors=F)
  colnames(df_files)<-c("sample","file")
  
  cat("df_files_0\n")
  cat(str(df_files))
  cat("\n")
    
    
 
  
  
  ##### Reading loop -----
  
  file_array<-df_files$file
  
  
  list_results<-list()
  

  
  DEBUG<-1
  
  for(i in 1:length(file_array))
  {
    file_array_sel<-file_array[i]
    
    cat("------------------------------------------------>\t")
    cat(sprintf(as.character(file_array_sel)))
    cat("\t")
    
    df_files_sel<-df_files[which(df_files$file == file_array_sel),]
    
    sample_sel<-df_files_sel$sample
    
    cat(sprintf(as.character(sample_sel)))
    cat("\n")
    
    if(DEBUG == 1)
    {
      cat("df_files_sel_0\n")
      cat(str(df_files_sel))
      cat("\n")
     
      
    }

    if(dim(df_files_sel)[1] >0)
    {

      ### Set directory
      
      setwd(indir)
      
      SIZE_gate<-file.info(file_array_sel)$size
      
      if(DEBUG == 1)
      {
        cat("SIZE_gate_0\n")
        cat(str(SIZE_gate))
        cat("\n")
      }
      
      if(SIZE_gate> 0)
      {
        LINE_gate<-length(readLines(file_array_sel))
        
        if(DEBUG == 1)
        {
          cat("LINE_gate_0\n")
          cat(str(LINE_gate))
          cat("\n")
        }
        
        if(LINE_gate> 0)
        {
          results<-as.data.frame(fread(file=file_array_sel, sep=" ", header = F, skip=1, fill=TRUE), stringsAsFactors = F)
          
          colnames(results)<-c("unknown","CellBC", "MolBC", "GFPbc")
          
          results$sample<-sample_sel
          
          results$GFPbc[which(results$GFPbc == "")]<-'No_GFPbcs'
          
          if(DEBUG == 1)
          {
            cat("results_0\n")
            cat(str(results))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(results$GFPbc))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(results$GFPbc)))))
            cat("\n")
          }
          
          # results$UMI_GFP = paste(results$MolBC, results$GFPbc, sep="-")
          # results$Cell_GFP = paste(results$CellBC, results$GFPbc, sep="-")
          results$UMI_GFP_Cell = paste(results$MolBC, results$GFPbc,results$CellBC, sep="-")
          
          if(DEBUG == 1)
          {
            cat("results_1\n")
            cat(str(results))
            cat("\n")
          }
          
          results<-results[,-which(colnames(results) == 'unknown')]
          
          if(DEBUG == 1)
          {
            cat("results_2\n")
            cat(str(results))
            cat("\n")
          }
          
          
          ## Remove multiple reads of the same UMI, Cell and GFP barcode (UMI duplicates)
          
          results_no_duplicated_UMIS <- results[!duplicated(results$UMI_GFP_Cell),]
          
          if(DEBUG == 1)
          {
            cat("results_no_duplicated_UMIS_0\n")
            cat(str(results_no_duplicated_UMIS))
            cat("\n")
          }
          
          ## find how many unique UMIS support each cell GFP combination ----
          
          results_no_duplicated_UMIS.dt<-data.table(results_no_duplicated_UMIS, key=c("sample","CellBC", "GFPbc"))
          
          
          
          results_no_duplicated_UMIS_collapsed_UMIS<-as.data.frame(results_no_duplicated_UMIS.dt[,.(Number_of_UMIS=.N), by=key(results_no_duplicated_UMIS.dt)], stringsAsFactors=F)
          
          
          if(DEBUG == 1)
          {
            cat("results_no_duplicated_UMIS_collapsed_UMIS_0\n")
            cat(str(results_no_duplicated_UMIS_collapsed_UMIS))
            cat("\n")
          }
          
          ## filter by Threshold_UMIS_per_cell ------------------
          
          
          results_no_duplicated_UMIS_collapsed_UMIS_PASS_filter_reads<-results_no_duplicated_UMIS_collapsed_UMIS[which(results_no_duplicated_UMIS_collapsed_UMIS$Number_of_UMIS >= Threshold_UMIS_per_cell),]
          
          cat("results_no_duplicated_UMIS_collapsed_UMIS_PASS_filter_reads_0\n")
          cat(str(results_no_duplicated_UMIS_collapsed_UMIS_PASS_filter_reads))
          cat("\n")
          cat(sprintf(as.character(names(summary(results_no_duplicated_UMIS_collapsed_UMIS_PASS_filter_reads$Number_of_UMIS)))))
          cat("\n")
          cat(sprintf(as.character(summary(results_no_duplicated_UMIS_collapsed_UMIS_PASS_filter_reads$Number_of_UMIS))))
          cat("\n")
          cat(sprintf(as.character(names(summary(as.factor(results_no_duplicated_UMIS_collapsed_UMIS_PASS_filter_reads$GFPbc))))))
          cat("\n")
          cat(sprintf(as.character(summary(as.factor(results_no_duplicated_UMIS_collapsed_UMIS_PASS_filter_reads$GFPbc)))))
          cat("\n")
          
          
          ## find cases where more than one genotype points to the same cell ----
          
          results_no_duplicated_UMIS_collapsed_UMIS_PASS_filter_reads.dt<-data.table(results_no_duplicated_UMIS_collapsed_UMIS_PASS_filter_reads, key=c("sample","CellBC"))
          
          
          
          Final_table<-as.data.frame(results_no_duplicated_UMIS_collapsed_UMIS_PASS_filter_reads.dt[,.(GFPbc_string=paste(GFPbc, collapse = ';'),
                                                                                                                   GFPbc_attribution=length(unlist(strsplit(paste(GFPbc, collapse = ';'), split=';'))),
                                                                                                                   Number_of_UMIS_string=paste(Number_of_UMIS, collapse = ';')), 
                                                                                                                by=key(results_no_duplicated_UMIS_collapsed_UMIS_PASS_filter_reads.dt)], stringsAsFactors=F)
          
          
          if(DEBUG == 1)
          {
            cat("Final_table_0\n")
            cat(str(Final_table))
            cat("\n")
            cat(sprintf(as.character(names(summary(Final_table$GFPbc_attribution)))))
            cat("\n")
            cat(sprintf(as.character(summary(Final_table$GFPbc_attribution))))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(Final_table$GFPbc_string))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(Final_table$GFPbc_string)))))
            cat("\n")
          }
          
          Final_table_uniquely_genotyped<-Final_table[which(Final_table$GFPbc_attribution == Threshold_attributed_genotypes),]
          

          Final_table_uniquely_genotyped<-Final_table_uniquely_genotyped[which(Final_table$GFPbc_attribution == Threshold_attributed_genotypes),]
          
          
          cat("Final_table_uniquely_genotyped_0\n")
          cat(str(Final_table_uniquely_genotyped))
          cat("\n")
          cat(str(unique(Final_table_uniquely_genotyped$CellBC)))
          cat("\n")
          cat(sprintf(as.character(names(summary(as.factor(Final_table_uniquely_genotyped$GFPbc_string))))))
          cat("\n")
          cat(sprintf(as.character(summary(as.factor(Final_table_uniquely_genotyped$GFPbc_string)))))
          cat("\n")
          

          list_results[[i]]<-Final_table_uniquely_genotyped
          
          
          
        }#LINE_gate> 0
      }# SIZE_gate> 0
    }# dim(df_files_sel)[1] >0
  }# i in 1:length(file_array)
  
  
  if(length(list_results) >0)
  {
    Results = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
    cat("Results_0\n")
    cat(str(Results))
    cat("\n")
    cat(str(unique(Results$sample)))
    cat("\n")
   
    
   
    
    colnames(Results)[which(colnames(Results) == 'Number_of_UMIS_string')]<-'Number_of_UMIS'
    colnames(Results)[which(colnames(Results) == 'GFPbc_string')]<-'GFPbc'
    
    Results$Number_of_UMIS<-as.integer(Results$Number_of_UMIS)
    
    Results$sample<-factor(Results$sample,
                                              levels=samples,
                                              ordered=T)
    Results$GFPbc<-factor(Results$GFPbc,
                                             levels=genotypes_string,
                                             ordered=T)
    
    Results<-Results[order(Results$sample,Results$GFPbc),]
    
    
    cat("Results_1\n")
    cat(str(Results))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results$Number_of_UMIS)))))
    cat("\n")
    cat(sprintf(as.character(summary(Results$Number_of_UMIS))))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results$GFPbc)))))
    cat("\n")
    cat(sprintf(as.character(summary(Results$GFPbc))))
    cat("\n")
    
    
    Results<-Results[!is.na(Results$GFPbc),]
    
    
    cat("Results_2\n")
    cat(str(Results))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results$Number_of_UMIS)))))
    cat("\n")
    cat(sprintf(as.character(summary(Results$Number_of_UMIS))))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results$GFPbc)))))
    cat("\n")
    cat(sprintf(as.character(summary(Results$GFPbc))))
    cat("\n")
    
   
    ### save ----
    
    setwd(out)
    
    write.table(Results, file='Uniquely_genotyped_and_filtered_larry_barcodes_assignments.csv', sep=",",quote=F,row.names=F)
    
    saveRDS(Results, file='Uniquely_genotyped_and_filtered_larry_barcodes_assignments.rds')
   
    
    
  }# length(list_results) >0
}

graph_function = function(option_list)
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
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  path_graphs<-paste(out,'graphs','/',sep='')
  
  if(file.exists(path_graphs)){
    
  }else{
    dir.create(file.path(path_graphs))
  }#path_graphs
  
  
  #### READ and transform Threshold_UMIS_per_cell ----
  
  Threshold_UMIS_per_cell = opt$Threshold_UMIS_per_cell
  
  cat("Threshold_UMIS_per_cell_\n")
  cat(sprintf(as.character(Threshold_UMIS_per_cell)))
  cat("\n")
  
  ### read results unfiltered ----
  
  setwd(out)
  
  Results_uniquely_genotyped<-readRDS(file='Uniquely_genotyped_and_filtered_larry_barcodes_assignments.rds')
  
  
  cat("Results_uniquely_genotyped_0\n")
  cat(str(Results_uniquely_genotyped))
  cat("\n")
  cat(sprintf(as.character(names(summary(Results_uniquely_genotyped$Number_of_UMIS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Results_uniquely_genotyped$Number_of_UMIS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Results_uniquely_genotyped$GFPbc)))))
  cat("\n")
  cat(sprintf(as.character(summary(Results_uniquely_genotyped$GFPbc))))
  cat("\n")
  
  #### LOOP to print the different files at increasing thresholds ----
  
  Threshold_array<-seq(Threshold_UMIS_per_cell,Threshold_UMIS_per_cell+2, by=1)
  
  cat("Threshold_array_0\n")
  cat(str(Threshold_array))
  cat("\n")
  
  DEBUG<-1
  
  for(i in 1:length(Threshold_array)){
    
    Threshold_array_sel<-Threshold_array[i]
    
    cat("Threshold------------------------------------------------------------------->\t")
    cat(sprintf(as.character(Threshold_array_sel)))
    cat("\n")
    
    Results_uniquely_genotyped_sel<-Results_uniquely_genotyped[which(Results_uniquely_genotyped$Number_of_UMIS >= Threshold_array_sel),]
    
    if(DEBUG == 1){
      
      cat("Results_uniquely_genotyped_sel_0\n")
      cat(str(Results_uniquely_genotyped_sel))
      cat("\n")
    }
    
    Results_uniquely_genotyped_sel.dt<-data.table(Results_uniquely_genotyped_sel, key=c("sample","GFPbc"))
    
    
    Freq_table<-as.data.frame(Results_uniquely_genotyped_sel.dt[,.(Freq=.N),by=key(Results_uniquely_genotyped_sel.dt)], stringsAsFactors=F)
    
    if(DEBUG == 1){
      
      cat("Freq_table_0\n")
      cat(str(Freq_table))
      cat("\n")
    }
    
    Results_uniquely_genotyped_sel.dt<-data.table(Results_uniquely_genotyped_sel, key=c("sample"))
    
    
    Freq_TOTAL<-as.data.frame(Results_uniquely_genotyped_sel.dt[,.(TOTAL=.N),by=key(Results_uniquely_genotyped_sel.dt)], stringsAsFactors=F)
    
    if(DEBUG == 1){
      
      cat("Freq_TOTAL_0\n")
      cat(str(Freq_TOTAL))
      cat("\n")
    }
    
    Results_uniquely_genotyped_sel.dt<-data.table(Results_uniquely_genotyped_sel, key=c("GFPbc"))
    
    
    Freq_GFPbc<-as.data.frame(Results_uniquely_genotyped_sel.dt[,.(TOTAL=.N),by=key(Results_uniquely_genotyped_sel.dt)], stringsAsFactors=F)
    
    if(DEBUG == 1){
      
      cat("Freq_GFPbc_0\n")
      cat(str(Freq_GFPbc))
      cat("\n")
    }
    
    Freq_table<-merge(Freq_table,
                      Freq_TOTAL,
                      by="sample")
    
    Freq_table$Perc<-round(100*(Freq_table$Freq/Freq_table$TOTAL),2)
    
    if(DEBUG == 1){
      
      cat("Freq_table_1\n")
      cat(str(Freq_table))
      cat("\n")
    }
    
    
    #### Stacked Graph ---------------
    
    breaks.Rank<-(seq(0,100,by=25))
    labels.Rank<-as.character(breaks.Rank)
    
    if(DEBUG == 1){
      cat("-------------------------------------->\t")
      cat(sprintf(as.character(labels.Rank)))
      cat("\n")
    }
    
    fill_colours<-c(brewer.pal(9, "Greens")[c(5,6,7)],brewer.pal(9, "Reds")[c(5,6,7)],brewer.pal(9, "Purples")[c(5,6,7)],brewer.pal(9, "Blues")[c(4,5,6)],'gray','black')
     
    
    if(DEBUG == 1){
      cat("--------------------fill_colours------------------>\t")
      cat(str(fill_colours))
      cat("\n")
      cat(sprintf(as.character(fill_colours)))
      cat("\n")
    }
    
    
    stacked_barplot<-Freq_table %>%
      mutate(myaxis = paste0(sample,"\n", "n=", TOTAL), drop=F) %>%
      mutate(myaxis=fct_reorder(myaxis,as.numeric(sample)), drop=F) %>%
      ggplot(aes(x=myaxis, y=Perc, fill=GFPbc)) +
      geom_bar(stat="identity",colour='white')+
      scale_y_continuous(name=paste("Percentage of total cells",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                         limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
      scale_x_discrete(name=NULL, drop=F)+
      scale_fill_manual(values=fill_colours,
                        drop=F,
                     name="GFPbc", breaks=Freq_GFPbc$GFPbc,
                     labels=paste(Freq_GFPbc$GFPbc,
                                  Freq_GFPbc$TOTAL, sep =' n= '))+
      ggtitle(paste("Min.",Threshold_array_sel,"different UMIs per GFPbc and cell ", sep=' '))+
      theme_classic()+
      theme(plot.title=element_text(color="black", family="sans"),
            axis.title.y=element_text(color="black", family="sans"),
            axis.title.x=element_blank(),
            axis.text.y=element_text(angle=0,color="black", family="sans"),
            axis.text.x=element_text(angle=45,vjust=1,hjust=1,color="black", family="sans"))+
      theme(legend.title = element_text(color="black", family="sans"),
            legend.text = element_text(color="black", family="sans"),
            legend.key.size = unit(0.35, 'cm'), #change legend key size
            legend.key.height = unit(0.35, 'cm'), #change legend key height
            legend.key.width = unit(0.35, 'cm'), #change legend key width
            legend.position="right")+
      ggeasy::easy_center_title()
      
      
    setwd(path_graphs)
    
    saveRDS(file=paste("Graph_","Threshold_",Threshold_array_sel,'.rds', sep=''), stacked_barplot)
    
    svgname<-paste(paste("test",Threshold_array_sel, sep='_'),".svg",sep='')
    svglite(svgname, width = 5, height = 3)
    print(stacked_barplot)
    dev.off()
  }#i in 1:length(Threshold_array)
  
}

patch_work_function = function(option_list)
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
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  path_graphs<-paste(out,'graphs','/',sep='')
  
  if(file.exists(path_graphs)){
    
  }else{
    dir.create(file.path(path_graphs))
  }#path_graphs
  
 
  
  ############## Retrieve graphs  -----------------
  
  
  file_list <- list.files(path=path_graphs, include.dirs = FALSE)
  
  
  cat("file_list\n")
  cat(str(file_list))
  cat("\n")
  
  
  indexes_sel <- grep("\\.rds",file_list)
  
  file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
  colnames(file_list_sel)<-"file"
  
  cat("file_list_sel_0\n")
  cat(str(file_list_sel))
  cat("\n")
  
  file_list_sel$filter<-gsub("Graph_Threshold_","",file_list_sel$file)
  file_list_sel$filter<-gsub("\\.rds","",file_list_sel$filter)
  
  
  cat("file_list_sel_1\n")
  cat(str(file_list_sel))
  cat("\n")  
  
  
  
  ############# LOOP --------------------------
  
  List_RESULTS<-list()
  
  
  # Results_DEF<-data.frame()
  
  DEBUG<-0
  
  for(i in 1:dim(file_list_sel)[1]){
    
    setwd(path_graphs)
    
    read_file_sel<-file_list_sel$file[i]
    filter_sel<-file_list_sel$filter[i]
    
    
    cat("-------------------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(read_file_sel)))
    cat("\t")
    cat(sprintf(as.character(filter_sel)))
    cat("\n")
    
    graph<-readRDS(file=read_file_sel)
    
   
    
    
    List_RESULTS[[i]]<-graph
    
  }#i in 1:dim(file_list_sel)[1]
  
  
  # cat("List_RESULTS_0\n")
  # cat(str(List_RESULTS))
  # cat("\n")
  
  
  ########################## SAVE PLOTS -------------------------------------------------
  
  setwd(path_graphs)
  
  p <- ((List_RESULTS[[1]] / List_RESULTS[[2]] / List_RESULTS[[3]])) +
    plot_annotation(tag_levels = "a") +
    plot_layout( guides = "collect") & theme(text = element_text(size = 8, family = "Arial"))
  
  ggsave(p, filename = "plot.pdf", device=cairo_pdf, height=20, width=16, unit="cm", path=NULL, scale=1)
  ggsave(p, filename = "plot.svg", device=svg, height=20, width=16, unit="cm", path=NULL, scale=1)
  
  
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
    make_option(c("--samples"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--genotypes_string"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--indir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_attributed_genotypes"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_UMIS_per_cell"), type="numeric", default=NULL, 
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
  
  read_all_barcoded_samples(opt)
  graph_function(opt)
  patch_work_function(opt)


  
}


###########################################################################

system.time( main() )
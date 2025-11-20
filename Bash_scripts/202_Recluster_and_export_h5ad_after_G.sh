#!/bin/bash

eval "$(conda shell.bash hook)"
 
  
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
Pythonscripts_path=$(echo "/home/manuel.tardaguila/Scripts/Python/")

MASTER_ROUTE=$1
analysis=$2


##########################################################################

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

mkdir -p $output_dir

Log_files=$(echo "$output_dir""/""Log_files/")
 
mkdir -p $Log_files
 
### Recluster at 

res_param=$(echo '0.5')

type=$(echo "Recluster_at_""$res_param")
outfile_Recluster=$(echo "$Log_files""outfile_11_""$type"".log")
touch $outfile_Recluster
echo -n "" > $outfile_Recluster
name_Recluster=$(echo "$type""_job")


Rscript_Recluster=$(echo "$Rscripts_path""506_Clustering_of_merged_samples_after_G.R")

db_filt_clustered_QCed_cell_annotated_MJ_only_genotyped=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_QCed_cell_annotated_MJ_only_genotyped.rds")
mem=$(echo "8192")
processors=$(echo "15")
total_memory=$(( mem * processors ))

echo "$processors"
echo "$total_memory"

# 15 8192
 
myjobid_Recluster=$(sbatch --job-name $name_Recluster --output=$outfile_Recluster --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/multiome_QC_DEF/ Rscript $Rscript_Recluster --db_filt_clustered_QCed_cell_annotated_MJ_only_genotyped $db_filt_clustered_QCed_cell_annotated_MJ_only_genotyped --res_param $res_param --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_Recluster=$(sbatch --dependency=afterany:$myjobid_Recluster --open-mode=append --output=$outfile_Recluster --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Recluster >> $outfile_Recluster")


### MACS2_calling_by_annotation 

type=$(echo "MACS2_calling_by_annotation")
outfile_MACS2_calling_by_annotation=$(echo "$Log_files""outfile_12_""$type"".log")
touch $outfile_MACS2_calling_by_annotation
echo -n "" > $outfile_MACS2_calling_by_annotation
name_MACS2_calling_by_annotation=$(echo "$type""_job")


Rscript_MACS2_calling_by_annotation=$(echo "$Rscripts_path""507_MACS2_recall_peaks_from_cell_type_NEW_annotation.R")


db_filt_clustered_QCed_cell_annotated_MJ_only_genotyped_reclustered=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_QCed_cell_annotated_MJ_only_genotyped_reclustered.rds")
frag_file=$(echo "$output_dir""merged.atac_fragments.tsv.gz")


mem=$(echo "8192")
processors=$(echo "30")
total_memory=$(( mem * processors ))

echo "$processors"
echo "$total_memory"

# 15 8192



myjobid_MACS2_calling_by_annotation=$(sbatch --dependency=afterany:$myjobid_Recluster --job-name $name_MACS2_calling_by_annotation --output=$outfile_MACS2_calling_by_annotation --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/multiome_QC_DEF/ Rscript $Rscript_MACS2_calling_by_annotation --db_filt_clustered_QCed_cell_annotated_MJ_only_genotyped_reclustered $db_filt_clustered_QCed_cell_annotated_MJ_only_genotyped_reclustered --frag_file $frag_file --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_MACS2_calling_by_annotation=$(sbatch --dependency=afterany:$myjobid_MACS2_calling_by_annotation --open-mode=append --output=$outfile_MACS2_calling_by_annotation --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_MACS2_calling_by_annotation >> $outfile_MACS2_calling_by_annotation")

### Convert to h5ad ----

type=$(echo "convert_to_h5ad")
outfile_convert_to_h5ad=$(echo "$Log_files""outfile_13_""$type"".log")
touch $outfile_convert_to_h5ad
echo -n "" > $outfile_convert_to_h5ad
name_convert_to_h5ad=$(echo "$type""_job")


Pythonscript_convert_to_h5ad=$(echo "$Pythonscripts_path""4_export_to_h5ad_ATAC_and_RNA_for_SIMBA_v2.py")


metadata_file=$(echo "$output_dir""final_cell_metadata.rds")
rna_matrix_file=$(echo "$output_dir""final_rna_sct_matrix.rds")
atac_matrix_file=$(echo "$output_dir""final_atac_matrix.rds")
output_name=$(echo "merged_clusters_final")

mem=$(echo "4096")
processors=$(echo "4")
total_memory=$(( mem * processors ))

echo "$processors"
echo "$total_memory"



myjobid_convert_to_h5ad=$(sbatch --dependency=afterany:$myjobid_MACS2_calling_by_annotation --job-name $name_convert_to_h5ad --output=$outfile_convert_to_h5ad --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/h5ad_exporter/ python $Pythonscript_convert_to_h5ad --rna-matrix-file $rna_matrix_file --atac-matrix-file $atac_matrix_file --metadata-file $metadata_file --base-output-name $output_name --cores $processors --memory $total_memory --output-dir $output_dir")
myjobid_seff_convert_to_h5ad=$(sbatch --dependency=afterany:$myjobid_convert_to_h5ad --open-mode=append --output=$outfile_convert_to_h5ad --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_convert_to_h5ad >> $outfile_convert_to_h5ad")


### Linking_peaks_to_genes 

type=$(echo "Linking_peaks_to_genes")
outfile_Linking_peaks_to_genes=$(echo "$Log_files""outfile_14_""$type"".log")
touch $outfile_Linking_peaks_to_genes
echo -n "" > $outfile_Linking_peaks_to_genes
name_Linking_peaks_to_genes=$(echo "$type""_job")


Rscript_Linking_peaks_to_genes=$(echo "$Rscripts_path""479_LinkPeaks_function_v2.R")


SeuratObject=$(echo "$output_dir""merged_clusters_final.rds")
DE_genes_1=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/Downstream_analysis/DE_per_identity/DE_results_Diff_MK.rds")
DE_genes_2=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/Downstream_analysis/DE_per_identity/DE_results_Diff_lymph.rds")


mem=$(echo "8192")
processors=$(echo "15")
total_memory=$(( mem * processors ))

echo "$processors"
echo "$total_memory"

# 15 8192



myjobid_Linking_peaks_to_genes=$(sbatch --dependency=afterany:$myjobid_MACS2_calling_by_annotation --job-name $name_Linking_peaks_to_genes --output=$outfile_Linking_peaks_to_genes --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/multiome_QC_DEF/ Rscript $Rscript_Linking_peaks_to_genes --SeuratObject $SeuratObject --DE_genes_1 $DE_genes_1 --DE_genes_2 $DE_genes_2 --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_Linking_peaks_to_genes=$(sbatch --dependency=afterany:$myjobid_Linking_peaks_to_genes --open-mode=append --output=$outfile_Linking_peaks_to_genes --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Linking_peaks_to_genes >> $outfile_Linking_peaks_to_genes")


### create_chromatin_assay 

type=$(echo "create_chromatin_assay")
outfile_create_chromatin_assay=$(echo "$Log_files""outfile_15_""$type"".log")
touch $outfile_create_chromatin_assay
echo -n "" > $outfile_create_chromatin_assay
name_create_chromatin_assay=$(echo "$type""_job")


Rscript_create_chromatin_assay=$(echo "$Rscripts_path""509_Linked_Peaks_to_chromatin_assay_in_Seurat_object.R")


SeuratObject=$(echo "$output_dir""merged_clusters_final.rds")
Linked_peak_to_DE_genes_ATAC_by_Construction_annotation=$(echo "$output_dir""Linked_peak_to_selected_genes_ATAC_by_Construction_annotation.rds")


mem=$(echo "8192")
processors=$(echo "10")
total_memory=$(( mem * processors ))

echo "$processors"
echo "$total_memory"



myjobid_create_chromatin_assay=$(sbatch --dependency=afterany:$myjobid_Linking_peaks_to_genes --job-name $name_create_chromatin_assay --output=$outfile_create_chromatin_assay --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/multiome_QC_DEF/ Rscript $Rscript_create_chromatin_assay --SeuratObject $SeuratObject --Linked_peak_to_DE_genes_ATAC_by_Construction_annotation $Linked_peak_to_DE_genes_ATAC_by_Construction_annotation --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_create_chromatin_assay=$(sbatch --dependency=afterany:$myjobid_create_chromatin_assay --open-mode=append --output=$outfile_create_chromatin_assay --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_create_chromatin_assay >> $outfile_create_chromatin_assay")

 

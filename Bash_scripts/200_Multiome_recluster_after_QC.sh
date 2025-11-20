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
outfile_Recluster=$(echo "$Log_files""outfile_8_""$type"".log")
touch $outfile_Recluster
echo -n "" > $outfile_Recluster
name_Recluster=$(echo "$type""_job")


Rscript_Recluster=$(echo "$Rscripts_path""504_Clustering_of_merged_samples_after_QC.R")

db_filt_clustered_QCed=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_QCed.rds")
mem=$(echo "8192")
processors=$(echo "20")
total_memory=$(( mem * processors ))

echo "$processors"
echo "$total_memory"

# 15 8192
 
myjobid_Recluster=$(sbatch --job-name $name_Recluster --output=$outfile_Recluster --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/multiome_QC_DEF/ Rscript $Rscript_Recluster --db_filt_clustered_QCed $db_filt_clustered_QCed --res_param $res_param --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_Recluster=$(sbatch --dependency=afterany:$myjobid_Recluster --open-mode=append --output=$outfile_Recluster --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Recluster >> $outfile_Recluster")



### Conver to h5ad ----

type=$(echo "convert_to_h5ad")
outfile_convert_to_h5ad=$(echo "$Log_files""outfile_9_""$type"".log")
touch $outfile_convert_to_h5ad
echo -n "" > $outfile_convert_to_h5ad
name_convert_to_h5ad=$(echo "$type""_job")


Pythonscript_convert_to_h5ad=$(echo "$Pythonscripts_path""3_export_to_h5ad_v2.py")


metadata_file=$(echo "$output_dir""final_cell_metadata.rds")
matrix_file=$(echo "$output_dir""final_sct_matrix.rds")
output_name=$(echo "merged_unprocessed_db_filt_clustered_QCed_reclustered.h5ad")

mem=$(echo "4096")
processors=$(echo "4")
total_memory=$(( mem * processors ))

echo "$processors"
echo "$total_memory"


myjobid_convert_to_h5ad=$(sbatch --dependency=afterany:$myjobid_Recluster --job-name $name_convert_to_h5ad --output=$outfile_convert_to_h5ad --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/h5ad_exporter/ python $Pythonscript_convert_to_h5ad --matrix-file $matrix_file --metadata-file $metadata_file --output-name $output_name --cores $processors --memory $total_memory --output-dir $output_dir")
myjobid_seff_convert_to_h5ad=$(sbatch --dependency=afterany:$myjobid_convert_to_h5ad --open-mode=append --output=$outfile_convert_to_h5ad --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_convert_to_h5ad >> $outfile_convert_to_h5ad")




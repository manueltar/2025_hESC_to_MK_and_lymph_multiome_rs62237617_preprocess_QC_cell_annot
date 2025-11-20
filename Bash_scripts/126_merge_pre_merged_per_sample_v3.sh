#!/bin/bash
 
eval "$(conda shell.bash hook)"
   
   
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2
frag_file=$3

##########################################################################

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 

sample_array=$(echo 'MCO_01326,MCO_01327,MCO_01328,MCO_01329,MCO_01330,MCO_01331,MCO_01332,MCO_01333')
 
### Merge_pre_merged_per_sample

type=$(echo "$sample_array_sel""_""Merge_pre_merged_per_sample")
outfile_Merge_pre_merged_per_sample=$(echo "$Log_files""outfile_6_""$type"".log")
touch $outfile_Merge_pre_merged_per_sample
echo -n "" > $outfile_Merge_pre_merged_per_sample
name_Merge_pre_merged_per_sample=$(echo "$type""_job")


Rscript_Merge_pre_merged_per_sample=$(echo "$Rscripts_path""406_Merge_samples_vFINAL.R")

mem=$(echo "8184")
processors=$(echo "30")
total_memory=$(( mem * processors ))


#### 30 8192

path_MASTER_ROUTE=$(echo "$output_dir")


myjobid_Merge_pre_merged_per_sample=$(sbatch --job-name $name_Merge_pre_merged_per_sample --output=$outfile_Merge_pre_merged_per_sample --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/multiome_QC_DEF/ Rscript $Rscript_Merge_pre_merged_per_sample --sample_array $sample_array --frag_file $frag_file --type $type --processors $processors --total_memory $total_memory --path_MASTER_ROUTE $path_MASTER_ROUTE")
myjobid_seff_Merge_pre_merged_per_sample=$(sbatch --dependency=afterany:$myjobid_Merge_pre_merged_per_sample --open-mode=append --output=$outfile_Merge_pre_merged_per_sample --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Merge_pre_merged_per_sample >> $outfile_Merge_pre_merged_per_sampley")

##########################################################################################################################

### cluster_merged_object

type=$(echo "cluster_merged_object")
outfile_cluster_merged_object=$(echo "$Log_files""outfile_7_""$type"".log")
touch $outfile_cluster_merged_object
echo -n "" > $outfile_cluster_merged_object
name_cluster_merged_object=$(echo "$type""_job")


Rscript_cluster_merged_object=$(echo "$Rscripts_path""408_Clustering_of_merged_samples_v2.R")

filtered_db_object=$(echo "$output_dir""merged_unprocessed_db_filt.rds")
mem=$(echo "8184")
processors=$(echo "30")
total_memory=$(( mem * processors ))


# --dependency=afterany:$myjobid_Merge_pre_merged_per_sample
 
myjobid_cluster_merged_object=$(sbatch --dependency=afterany:$myjobid_Merge_pre_merged_per_sample --job-name $name_cluster_merged_object --output=$outfile_cluster_merged_object --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/multiome_QC_DEF/ Rscript $Rscript_cluster_merged_object --filtered_db_object $filtered_db_object --type $type --out $output_dir --processors $processors --total_memory $total_memory")
myjobid_seff_cluster_merged_object=$(sbatch --dependency=afterany:$myjobid_cluster_merged_object --open-mode=append --output=$outfile_cluster_merged_object --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_cluster_merged_object >> $outfile_cluster_merged_object")




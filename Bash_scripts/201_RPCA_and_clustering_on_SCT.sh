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



type=$(echo "RPCA_and_recluster_at_""$res_param")
outfile_RPCA_and_recluster=$(echo "$Log_files""outfile_16_""$type"".log")
touch $outfile_RPCA_and_recluster
echo -n "" > $outfile_RPCA_and_recluster
name_RPCA_and_recluster=$(echo "$type""_job")


Rscript_RPCA_and_recluster=$(echo "$Rscripts_path""505_rpca_integration_and_clustering_multiome.R")
SeuratObject=$(echo "$output_dir""merged_clusters_final_with_peak_subset.rds")

mem=$(echo "8192")
processors=$(echo "20")
total_memory=$(( mem * processors ))

echo "$processors"
echo "$total_memory"


 
myjobid_RPCA_and_recluster=$(sbatch --job-name $name_RPCA_and_recluster --output=$outfile_RPCA_and_recluster --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/multiome_QC_DEF/ Rscript $Rscript_RPCA_and_recluster --SeuratObject $SeuratObject --res_param $res_param --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_RPCA_and_recluster=$(sbatch --dependency=afterany:$myjobid_RPCA_and_recluster --open-mode=append --output=$outfile_RPCA_and_recluster --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_RPCA_and_recluster >> $outfile_RPCA_and_recluster")




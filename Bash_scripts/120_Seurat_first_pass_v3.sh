#!/bin/bash

eval "$(conda shell.bash hook)"
  

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

mkdir -p $output_dir

Log_files=$(echo "$output_dir""/""Log_files/")
 
mkdir -p $Log_files


sample_array=$(echo 'MCO_01326,MCO_01327,MCO_01328,MCO_01329,MCO_01330,MCO_01331,MCO_01332,MCO_01333')
sample_array=$(echo 'MCO_01327,MCO_01328,MCO_01329,MCO_01330,MCO_01331,MCO_01332,MCO_01333')

 
a=($(echo "$sample_array" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
do

    sample_array_sel=$i
    echo "$sample_array_sel"

 
    ### Seurat_first_pass

    type=$(echo "$sample_array_sel""_""Seurat_first_pass")
    outfile_Seurat_first_pass=$(echo "$Log_files""outfile_1_""$type"".log")
    touch $outfile_Seurat_first_pass
    echo -n "" > $outfile_Seurat_first_pass
    name_Seurat_first_pass=$(echo "$type""_job")


    Rscript_Seurat_first_pass=$(echo "$Rscripts_path""403_Seurat_first_pass_v3.R")

    sample_name=$sample_array_sel
    master_path=$MASTER_ROUTE
    rna_min_features=$(echo '500')
    atac_min_fragments=$(echo '1000')
    MITO_max=$(echo '10')
    mem=$(echo "8184")
    processors=$(echo "16")
    total_memory=$(( mem * processors ))
    
    myjobid_Seurat_first_pass=$(sbatch --job-name $name_Seurat_first_pass --output=$outfile_Seurat_first_pass --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/multiome_QC_DEF/ Rscript $Rscript_Seurat_first_pass --sample_name $sample_name --master_path $master_path --rna_min_features $rna_min_features --atac_min_fragments $atac_min_fragments --MITO_max $MITO_max  --processors $processors --total_memory $total_memory --type $type --out $output_dir")
    myjobid_seff_Seurat_first_pass=$(sbatch --dependency=afterany:$myjobid_Seurat_first_pass --open-mode=append --output=$outfile_Seurat_first_pass --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Seurat_first_pass >> $outfile_Seurat_first_pass")


done



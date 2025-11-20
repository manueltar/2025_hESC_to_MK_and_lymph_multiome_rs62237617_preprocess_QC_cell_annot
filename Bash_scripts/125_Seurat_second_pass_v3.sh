#!/bin/bash

eval "$(conda shell.bash hook)"

MASTER_ROUTE=$1


Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")


sample_array=$(echo 'MCO_01326,MCO_01327,MCO_01328,MCO_01329,MCO_01330,MCO_01331,MCO_01332,MCO_01333')


a=($(echo "$sample_array" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
do

    sample_array_sel=$i
    echo "$sample_array_sel"



     preliminary_filtered=$(echo "$MASTER_ROUTE""processing_outputs""/""$sample_array_sel""/""intermediate""/""preliminary_filtered.rds")
     path_processing_outputs=$(echo "$MASTER_ROUTE""processing_outputs""/""$sample_array_sel""/")
     intermediate_dir=$(echo "$MASTER_ROUTE""processing_outputs""/""$sample_array_sel""/""intermediate""/")
     snATAC_dir=$(echo "$MASTER_ROUTE""processing_outputs""/""$sample_array_sel""/""snATAC_matrices""/")
     crange_dir=$(echo "$MASTER_ROUTE""$sample_array_sel""/""outs""/")
     premerge_dir=$(echo "$MASTER_ROUTE""processing_outputs""/""$sample_array_sel""/""pre_merge""/")
     fragfile=$(echo "$MASTER_ROUTE""$sample_array_sel""/""outs""/""atac_fragments.tsv.gz")

     output_dir=$(echo "$MASTER_ROUTE""processing_outputs""/")
     Log_files=$(echo "$output_dir""/""Log_files/")

   echo "$preliminary_filtered"
   echo "$path_processing_outputs"
   echo "$intermediate_dir"
   echo "$snATAC_dir"
   echo "$crange_dir"
   echo "$premerge_dir"
   echo "$fragfile"

   

   mkdir -p $premerge_dir

    ### Seurat_second_pass

    type=$(echo "$sample_array_sel""_""Seurat_second_pass")
    outfile_Seurat_second_pass=$(echo "$Log_files""outfile_5_""$type"".log")
    touch $outfile_Seurat_second_pass
    echo -n "" > $outfile_Seurat_second_pass
    name_Seurat_second_pass=$(echo "$type""_job")

 
    Rscript_Seurat_second_pass=$(echo "$Rscripts_path""405_Seurat_second_pass_v5.R")

    sample_name=$sample_array_sel
    mem=$(echo "8184")
    processors=$(echo "15")
    total_memory=$(( mem * processors ))

 

    myjobid_Seurat_second_pass=$(sbatch --job-name $name_Seurat_second_pass --output=$outfile_Seurat_second_pass --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/multiome_QC_DEF/ Rscript $Rscript_Seurat_second_pass --sample_name $sample_name --preliminary_filtered $preliminary_filtered --path_processing_outputs $path_processing_outputs --intermediate_dir $intermediate_dir --snATAC_dir $snATAC_dir --crange_dir $crange_dir --premerge_dir $premerge_dir --fragfile $fragfile --type $type --out $output_dir  --processors $processors --total_memory $total_memory")
    myjobid_seff_Seurat_second_pass=$(sbatch --dependency=afterany:$myjobid_Seurat_second_pass --open-mode=append --output=$outfile_Seurat_second_pass --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Seurat_second_pass >> $outfile_Seurat_second_pass")


done



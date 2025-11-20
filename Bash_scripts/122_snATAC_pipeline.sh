#!/bin/bash

eval "$(conda shell.bash hook)"
 

python_path=$(echo "/home/manuel.tardaguila/Scripts/Python/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 

module load samtools
module load bedtools
module load R/4.3.1

#conda activate /home/manuel.tardaguila/conda_envs/Manuel_ATAC
 

#sample_array=$(echo 'MCO_01326,MCO_01327,MCO_01328,MCO_01329')
#sample_array=$(echo 'MCO_01330,MCO_01331,MCO_01332,MCO_01333')
#sample_array=$(echo 'MCO_1278,MCO_1279,MCO_1280,MCO_1281')


sample_array=$(echo 'MCO_01326,MCO_01327,MCO_01328,MCO_01329,MCO_01330,MCO_01331,MCO_01332,MCO_01333')

a=($(echo "$sample_array" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
do

    sample_array_sel=$i
    echo "$sample_array_sel"
 

    ### snATAC_pipeline

    type=$(echo "$sample_array_sel""_""snATAC_pipeline")
    outfile_snATAC_pipeline=$(echo "$Log_files""outfile_3_""$type"".log")
    touch $outfile_snATAC_pipeline
    echo -n "" > $outfile_snATAC_pipeline
    name_snATAC_pipeline=$(echo "$type""_job")


    python_script_snATAC_pipeline=$(echo "$python_path""snATAC_10x_matrices_pipeline_v2.py")

    sample_name=$sample_array_sel
    input_bam=$(echo "$MASTER_ROUTE""$sample_array_sel""/""outs/""atac_possorted_bam.bam")
    output_dir_of_atac=$(echo "$output_dir""$sample_array_sel""/""snATAC_matrices""/")
    barcode=$(echo "$output_dir""$sample_array_sel""/""intermediate/""keep_barcodes_step1.txt")


    mkdir -p $output_dir_of_atac

    echo "$input_bam"

    echo "$output_dir_of_atac"

    echo "$barcode"

    memory=$(echo '12000')    
    processors=$(echo '30')
    total_memory=$(( (memory/1000) * processors ))
    memory_gb=$(( memory/1500 ))

    echo "$memory_gb"
    echo "$total_memory"


    cd $output_dir_of_atac

    myjobid_snATAC_pipeline=$(sbatch --job-name $name_snATAC_pipeline --output=$outfile_snATAC_pipeline --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$memory --parsable --wrap="conda run -p /home/manuel.tardaguila/conda_envs/Manuel_ATAC/ python3 $python_script_snATAC_pipeline -m $memory_gb -t $processors -n $name_snATAC_pipeline -b $input_bam -o $output_dir_of_atac --keep-bc $barcode")
    myjobid_seff_snATAC_pipeline=$(sbatch --dependency=afterany:$myjobid_snATAC_pipeline --open-mode=append --output=$outfile_snATAC_pipeline --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_snATAC_pipeline >> $outfile_snATAC_pipeline")



done

#conda deactivate

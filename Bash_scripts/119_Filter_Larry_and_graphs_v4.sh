#!/bin/bash
    
MASTER_ROUTE=$1
analysis=$2
indir=$3


Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0


bashrc_file=$(echo "/home/manuel.tardaguila/.bashrc")

source $bashrc_file
eval "$(conda shell.bash hook)"
  

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

mkdir -p $output_dir

Log_files=$(echo "$output_dir""/""Log_files/")

mkdir -p $Log_files
    
#### Larry_counts_and_filter #############################

 
type=$(echo "Larry_counts_and_filter")
outfile_Larry_counts_and_filter=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_Larry_counts_and_filter
echo -n "" > $outfile_Larry_counts_and_filter
name_Larry_counts_and_filter=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")
 
Rscript_Larry_counts_and_filter=$(echo "$Rscripts_path""402_Larry_counts_CHEK2_final_version_with_Paolas_correction.R")


samples=$(echo "MCO_01326,MCO_01327,MCO_01328,MCO_01329,MCO_01330,MCO_01331,MCO_01332,MCO_01333")


Threshold_attributed_genotypes=$(echo "1")
Threshold_UMIS_per_cell=$(echo "3")
genotypes_string=$(echo "chrGFP_WTA,chrGFP_WTB,chrGFP_WTC,chrGFP_rs1,chrGFP_rs2,chrGFP_rs3,chrGFP_R882H1,chrGFP_R882H2,chrGFP_R882H3,chrGFP_rs_R882H1,chrGFP_rs_R882H2,chrGFP_rs_R882H3,No_GFPbcs")

myjobid_Larry_counts_and_filter=$(sbatch --job-name=$name_Larry_counts_and_filter --output=$outfile_Larry_counts_and_filter --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Larry_counts_and_filter --indir $indir --samples $samples --Threshold_attributed_genotypes $Threshold_attributed_genotypes --Threshold_UMIS_per_cell $Threshold_UMIS_per_cell --genotypes_string $genotypes_string --type $type --out $output_dir")
myjobid_seff_Larry_counts_and_filter=$(sbatch --dependency=afterany:$myjobid_Larry_counts_and_filter --open-mode=append --output=$outfile_Larry_counts_and_filter --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Larry_counts_and_filter >> $outfile_Larry_counts_and_filter")




#!/bin/bash
#SBATCH --job-name=name_Run_cellranger_arc_count_MCO_01330
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.tardaguila@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=outfile_Run_cellranger_arc_count_%j.log
#SBATCH --mem=64G
#SBATCH --time=36:00:00
 
output_dir=$1
sample_sel=$2

config_file=$(echo "$output_dir""$sample_sel""_config_file.csv")
echo "$config_file"

path_for_reference_for_alingment=$(echo "/group/soranzo/paola.benaglio/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0")

echo "$path_for_reference_for_alingment"


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"


module load cellranger-arc/2.0.2

export PATH=/home/paola.benaglio/cellranger-arc-2.0.2:$PATH

cd $output_dir
cellranger-arc count --id=$sample_sel \
               --reference=$path_for_reference_for_alingment \
               --libraries=$config_file \
               --localcores=16 \
               --localmem=64 \
               --jobmode=local

	

echo "========================"
echo "Completed: $(date)"




eval "$(conda shell.bash hook)"



#### Run_cellranger_arc_count in SLURM (https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/advanced/cluster-mode)

#!/bin/bash
#SBATCH --job-name=name_Run_cellranger_count
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.tardaguila@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=outfile_Run_cellranger_count_%j.log
#SBATCH --mem=64G
#SBATCH --time=36:00:00
   
output_dir=$1
sample_sel=$2
fastq_path=$3

path_for_10X_transcriptome_reference=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/GFP_transgene_vCHEK2_and_DNMT3A_cellranger")

echo "$path_for_10X_transcriptome_reference"


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

module load cellranger/8.0.1

cd $output_dir

sample_dir=$("$output_dir""$sample_sel""/")

rm -rf $sample_dir


cellranger count --id=$sample_sel \
	       --sample=$sample_sel \
               --transcriptome=$path_for_10X_transcriptome_reference \
               --fastqs=$fastq_path \
	       --chemistry ARC-v1 \
               --localcores=16 \
               --localmem=47 \
	       --create-bam=true \
               --jobmode=local

	

echo "========================"
echo "Completed: $(date)"




eval "$(conda shell.bash hook)"



#### Run_cellranger_count in SLURM (https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/advanced/cluster-mode)

#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#
 
#SBATCH --job-name=Cellbender
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.tardaguila@fht.org
#SBATCH --partition=gpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=Cellbender_%j.log
#SBATCH --mem=16G
#SBATCH --gres=gpu:1



echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

eval "$(conda shell.bash hook)" 


module load singularity
module load cuda

MASTER_ROUTE=$1

echo $MASTER_ROUTE

samp=$2

echo $samp
 

mkdir -p $MASTER_ROUTE/'processing_outputs'/$samp

cd $MASTER_ROUTE/'processing_outputs'/$samp

pwd

input_file=$(echo "$MASTER_ROUTE/$samp/outs/raw_feature_bc_matrix.h5")

echo $input_file

# $MASTER_ROUTE/'processing_outputs'/$samp/cellbender_gex.h5

singularity exec --nv --cleanenv -B /localscratch -B /group/soranzo/manuel.tardaguila/ -B /scratch /ssu/gassu/singularity/cellbender_0.3.0.sif cellbender remove-background \
	    --cuda \
	    --checkpoint $MASTER_ROUTE/'processing_outputs'/$samp/ckpt.tar.gz \
	    --input $input_file \
	    --output cellbender_gex.h5 \
	    --exclude-feature-types Peaks

conda activate cellbender_0.2.2

ptrepack --complevel 5 cellbender_gex.h5:/matrix cellbender_gex_seurat.h5:/matrix

conda deactivate

echo "========================"
echo "Completed: $(date)"

 

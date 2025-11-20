#!/bin/bash
#SBATCH --job-name=Deconvolve_Larry
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.tardaguila@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=Deconvolve_Larry_%j.log
#SBATCH --mem=64G
#SBATCH --time=36:00:00

module load samtools
module load bwa-mem2

output_dir=$1
special_reference=$2
samp=$3


Larry_deconvolve_dir=$(echo "$output_dir""deconvolute_LARRY")

mkdir -p $Larry_deconvolve_dir

cd $Larry_deconvolve_dir

echo "$samp"

dir=$(echo "$output_dir""$samp""/""outs")

samtools view -b -f 4 $dir/gex_possorted_bam.bam -@ 32 > $samp.unmapped.bam
samtools sort -@ 32 -n -o $samp.unmapped.sorted.bam $samp.unmapped.bam


samtools fastq -@ 32 $samp.unmapped.sorted.bam -T CB,UB,xf > $samp.unmapped.fq
sed -i 's/\t/_/g' $samp.unmapped.fq


### Remap to modified reference with deletions 
bwa-mem2 mem -M -t 16 $special_reference $samp.unmapped.fq | samtools sort -o $samp.realigned.bam - 
samtools index $samp.realigned.bam

### Remove multimappers to keep only reads uniquely mapping to each GFP barcode

samtools view -q 5 -b $samp.realigned.bam > $samp.realigned.unique.bam 
samtools index $samp.realigned.unique.bam


for geno_sel in  chrGFP_WTA chrGFP_WTB chrGFP_WTC chrGFP_rs1 chrGFP_rs2 chrGFP_rs3 chrGFP_R882H1 chrGFP_R882H2 chrGFP_R882H3 chrGFP_rs_R882H1 chrGFP_rs_R882H2 chrGFP_rs_R882H3;

do

    echo Genotype $geno_sel

    samtools view $samp.realigned.unique.bam $geno_sel|cut -f1 | sed 's/_/ /g' | awk '{print $3}' | sort | uniq > $samp.$geno_sel.unique.bcs
    samtools view $samp.realigned.unique.bam $geno_sel|cut -f1 | sed 's/_/ /g' | awk -v var="$geno_sel" '{print $3,$4,var}' | sort | uniq > $samp.$geno_sel.unique_bcs_umi
    samtools view $samp.realigned.unique.bam $geno_sel|cut -f1 | sed 's/_/ /g' | awk -v var="$geno_sel" '{print $2,$3,$4,var}' | sort  > $geno_sel.$samp.bcs_umi
done

cat *.$samp.bcs_umi | sort  > $samp.all_geno_bc_umi


echo "========================"
echo "Completed: $(date)"




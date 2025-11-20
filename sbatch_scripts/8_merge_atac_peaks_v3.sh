#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================


#SBATCH --job-name=merge_frag
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.tardaguila@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=merge_frag_%j.log
#SBATCH --mem=96G

# 16 96 --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"


 
module load htslib-tools

output_dir=$1

cd "$output_dir"


for SAMPLE in MCO_01326 MCO_01327 MCO_01328 MCO_01329 MCO_01330 MCO_01331 MCO_01332 MCO_01333;do


    echo "Processing: $SAMPLE $(date)" >&2
   

   FILE=$(echo "$output_dir""$SAMPLE""/outs/atac_fragments.tsv.gz")

   if [[ ! -f "$FILE" ]]; then
        echo "WARNING: File not found for $SAMPLE at $FILE" >&2
        continue
    fi

   echo "Using file: $FILE" >&2

 zcat "$FILE" | awk -v SAMPLE="$SAMPLE"  'BEGIN{FS=OFS="\t"} {print $1,$2,$3,SAMPLE"_"$4,$5}'
done | awk 'NF >= 5 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/' | sort -k1,1 -k2,2n -S 72G | bgzip -c -@ 16 > $output_dir/processing_outputs/merged.atac_fragments.tsv.gz 

tabix -p bed $output_dir/processing_outputs/merged.atac_fragments.tsv.gz




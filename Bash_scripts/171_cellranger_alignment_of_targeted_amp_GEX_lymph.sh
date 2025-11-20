#!/bin/bash

eval "$(conda shell.bash hook)"
 

MASTER_ROUTE=$1
analysis=$2

path_fastq=$3

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")


mkdir -p $Log_files

cell_ranger_dir=$(echo "$output_dir""/""cellranger/")


mkdir -p $cell_ranger_dir

module load trimmomatic/0.39
module load seqtk/1.4
module load seqkit/2.8.2
module load cellranger/6.1.2

trimmomatic_bin=$(echo "/share/apps/spack/latest-2024-01-30/linux-centos8-skylake_avx512/gcc-10.4.0/trimmomatic-0.39-3t2b7ls4ug3du4jdzyl3ur5md2q5yb2n/bin/trimmomatic-0.39.jar")



sample_array=$(echo 'MCO_01331,MCO_01332,MCO_01333')
equivalence_array=$(echo 'MCO_01374_LIB_S109_L002,MCO_01375_LIB_S110_L002,MCO_01376_LIB_S111_L002')

a=($(echo "$sample_array" | tr "," '\n'))
b=($(echo "$equivalence_array" | tr "," '\n'))

declare -a arr

array_2_length=${#a[@]}

for (( i=0; i<${array_2_length}; i=i+1 ));
do

    sample_array_sel=${a[$i]}
    echo "$sample_array_sel"


    equivalence_sel=${b[$i]}

    echo "$equivalence_sel"

    cell_ranger_sample=$(echo $equivalence_sel|sed -r 's/MCO_[0-9]+_LIB_//g')

    echo "1 $cell_ranger_sample"

    cell_ranger_sample=$(echo $cell_ranger_sample|sed -r 's/_L.+$//g')

    echo "2 $cell_ranger_sample"

    cell_ranger_sample=$(echo "$sample_array_sel""_""$cell_ranger_sample""_""L001")

    echo "3 $cell_ranger_sample"


    r1_L1=$(echo "$path_fastq""$equivalence_sel""_R1_001.fastq.gz")


    r1_L1_CROPPED=$(echo "$output_dir""$sample_array_sel""_R1_CROPPED.fastq.gz")
    r1_L1_CROPPED_FILTERED=$(echo "$output_dir""$sample_array_sel""_R1_CROPPED_FILTERED.fastq.gz")
    r1_FILTERED_list=$(echo "$output_dir""$sample_array_sel""_list.txt")

    r2_L1=$(echo "$path_fastq""$equivalence_sel""_R2_001.fastq.gz")

    r2_L1_HEADCROPPED=$(echo "$output_dir""$sample_array_sel""_R2_HEADCROPPED.fastq.gz")
    r2_L1_HEADCROPPED_FILTERED=$(echo "$output_dir""$sample_array_sel""_R2_HEADCROPPED_FILTERED.fastq")
    r2_L1_HEADCROPPED_FILTERED_zip=$(echo "$output_dir""$sample_array_sel""_R2_HEADCROPPED_FILTERED.fastq.gz")


 
     ################################################ TRIM R1 to meet exactly the format of 10X

     type=$(echo "$sample_array_sel""_crop_R1")
     outfile_crop_R1=$(echo "$Log_files""outfile_2_""$type"".log")
     touch $outfile_crop_R1
     echo -n "" > $outfile_crop_R1
     name_crop_R1=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")


    myjobid_crop_R1=$(sbatch --job-name=$name_crop_R1 --output=$outfile_crop_R1 --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="java -jar $trimmomatic_bin SE -threads 4 $r1_L1 $r1_L1_CROPPED CROP:28")
    myjobid_seff_crop_R1=$(sbatch --dependency=afterany:$myjobid_crop_R1 --open-mode=append --output=$outfile_crop_R1 --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_crop_R1 >> $outfile_crop_R1")


    ################################################ Filter R1 to be exactly 28 

     type=$(echo "$sample_array_sel""_filter_R1")
     outfile_filter_R1=$(echo "$Log_files""outfile_3_""$type"".log")
     touch $outfile_filter_R1
     echo -n "" > $outfile_filter_R1
     name_filter_R1=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")




    myjobid_filter_R1=$(sbatch --dependency=afterany:$myjobid_crop_R1 --job-name=$name_filter_R1 --output=$outfile_filter_R1 --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="seqkit seq -m 28 $r1_L1_CROPPED | gzip > $r1_L1_CROPPED_FILTERED")
    myjobid_seff_filter_R1=$(sbatch --dependency=afterany:$myjobid_filter_R1 --open-mode=append --output=$outfile_filter_R1 --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_filter_R1 >> $outfile_filter_R1")

    #  ################################################ list R1 IDs filtered 

     type=$(echo "$sample_array_sel""_list_R1_IDs")
     outfile_list_R1_IDs=$(echo "$Log_files""outfile_4_""$type"".log")
     touch $outfile_list_R1_IDs
     echo -n "" > $outfile_list_R1_IDs
     name_list_R1_IDs=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")


     # --dependency=afterany:$myjobid_filter_R1

    myjobid_list_R1_IDs=$(sbatch --dependency=afterany:$myjobid_filter_R1 --job-name=$name_list_R1_IDs --output=$outfile_list_R1_IDs --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="seqkit fx2tab $r1_L1_CROPPED_FILTERED | cut -d \" \" -f1 > $r1_FILTERED_list")
    myjobid_seff_list_R1_IDs=$(sbatch --dependency=afterany:$myjobid_list_R1_IDs --open-mode=append --output=$outfile_list_R1_IDs --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_list_R1_IDs >> $outfile_list_R1_IDs")



    ##### R2 ####################################################################################################################################################################################################################################################################################################
    ##### R2 ####################################################################################################################################################################################################################################################################################################
    ##### R2 ####################################################################################################################################################################################################################################################################################################

     type=$(echo "$sample_array_sel""_trim_R2")
     outfile_crop_R2=$(echo "$Log_files""outfile_6_""$type"".log")
     touch $outfile_crop_R2
     echo -n "" > $outfile_crop_R2
     name_crop_R2=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")


 
     myjobid_crop_R2=$(sbatch --job-name=$name_crop_R2 --output=$outfile_crop_R2 --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="java -jar $trimmomatic_bin SE -threads 4 $r2_L1 $r2_L1_HEADCROPPED HEADCROP:12")
     myjobid_seff_crop_R2=$(sbatch --dependency=afterany:$myjobid_crop_R2 --open-mode=append --output=$outfile_crop_R2 --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_crop_R2 >> $outfile_crop_R2")


     ################################################ Filter R2

     type=$(echo "$sample_array_sel""_FILTER_R2")
     outfile_FILTER_R2=$(echo "$Log_files""outfile_7_""$type"".log")
     touch $outfile_FILTER_R2
     echo -n "" > $outfile_FILTER_R2
     name_FILTER_R2=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")

     # --dependency=afterany:$myjobid_crop_R2:$myjobid_list_R1_IDs


     myjobid_FILTER_R2=$(sbatch --dependency=afterany:$myjobid_crop_R2:$myjobid_list_R1_IDs --job-name=$name_FILTER_R2 --output=$outfile_FILTER_R2 --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="seqtk subseq $r2_L1_HEADCROPPED $r1_FILTERED_list | gzip > $r2_L1_HEADCROPPED_FILTERED_zip")
     myjobid_seff_FILTER_R2=$(sbatch --dependency=afterany:$myjobid_FILTER_R2 --open-mode=append --output=$outfile_FILTER_R2 --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_FILTER_R2 >> $outfile_FILTER_R2")

     ########################################### MV_TO_CR_R1 #############################################################################################################################################################################################################################
     
     type=$(echo "$sample_array_sel""_MV_TO_CR_R1")
     outfile_MV_TO_CR_R1=$(echo "$Log_files""outfile_8_""$type"".log")
     touch $outfile_MV_TO_CR_R1
     echo -n "" > $outfile_MV_TO_CR_R1
     name_MV_TO_CR_R1=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")



     R1_dest=$(echo "$cell_ranger_dir""$cell_ranger_sample""_R1_001.fastq.gz")

     #  --dependency=afterany:$myjobid_filter_R1

     myjobid_MV_TO_CR_R1=$(sbatch --dependency=afterany:$myjobid_filter_R1 --job-name=$name_MV_TO_CR_R1 --output=$outfile_MV_TO_CR_R1 --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="scp -r $r1_L1_CROPPED_FILTERED $R1_dest")
     myjobid_seff_MV_TO_CR_R1=$(sbatch --dependency=afterany:$myjobid_MV_TO_CR_R1 --open-mode=append --output=$outfile_MV_TO_CR_R1 --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_MV_TO_CR_R1 >> $outfile_MV_TO_CR_R1")


     ########################################### MV_TO_CR_R2 #############################################################################################################################################################################################################################
     
     type=$(echo "$sample_array_sel""_MV_TO_CR_R2")
     outfile_MV_TO_CR_R2=$(echo "$Log_files""outfile_9_""$type"".log")
     touch $outfile_MV_TO_CR_R2
     echo -n "" > $outfile_MV_TO_CR_R2
     name_MV_TO_CR_R2=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")



     R2_dest=$(echo "$cell_ranger_dir""$cell_ranger_sample""_R2_001.fastq.gz")

     myjobid_MV_TO_CR_R2=$(sbatch --dependency=afterany:$myjobid_FILTER_R2 --job-name=$name_MV_TO_CR_R2 --output=$outfile_MV_TO_CR_R2 --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="scp -r $r2_L1_HEADCROPPED_FILTERED_zip $R2_dest")
     myjobid_seff_MV_TO_CR_R2=$(sbatch --dependency=afterany:$myjobid_MV_TO_CR_R2 --open-mode=append --output=$outfile_MV_TO_CR_R2 --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_MV_TO_CR_R2 >> $outfile_MV_TO_CR_R2")



done



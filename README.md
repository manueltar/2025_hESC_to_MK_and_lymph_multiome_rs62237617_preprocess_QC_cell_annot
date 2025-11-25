You got it\! Here is the ASCII diagram of the pipeline flow for easy copying:

```
+-------------------------------------------------------------+
|                                                             |
| 🔬 PHASE I: Targeted GEX Alignment & Correction             |
|                                                             |
+-------------------------------------------------------------+
       |
       v
+------------------------------------+
|  STEP 1: Custom Genome Indexing    |  <-- GFP_transgene_vCHEK2...
|  (cellranger mkref)                |
+------------------------------------+
       |
       v
+------------------------------------+
|  STEP 2: Targeted Read Adaptation  |  <-- Wraper_scripts/127/146/171
|  (Format FastQ for CellRanger)     |
+------------------------------------+
       | (Parallel for MCO_01326-01333)
       v
+------------------------------------+
|  STEP 3: CellRanger GEX Counting   |  <-- sbatch 18_Cell_ranger_count
|  (Alignment to custom reference)   |
+------------------------------------+
       |
       v
+------------------------------------+
|  STEP 4: Background Correction     |  <-- sbatch 15_cellbender
|  (CellBender on Targeted GEX)      |
+------------------------------------+
       |
       |
       +=====================================================+
                                |
                                v
+-------------------------------------------------------------+
|                                                             |
| 🧬 PHASE II: Multiome Integration & Annotation               |
|                                                             |
+-------------------------------------------------------------+
       |
       v
+------------------------------------+
|  BLOCK 1: CellRanger ArcCount      |
|  (Raw Multiome GEX/ATAC matrices)  |
+------------------------------------+
       |
       v
+------------------------------------+
|  BLOCK 2: Pre-processing & Alignment | <-- Initial Seurat, CellBender, snATAC,
|  (QC, Peak Merge, Barcode Align)     |     LARRY filtering
+------------------------------------+
       |
       v
+------------------------------------+
|  BLOCK 3: Doublet Detection (Amulet)|
+------------------------------------+
       |
       v
+------------------------------------+
|  BLOCK 4 & 5: Seurat Re-pass & Merge |
|  (Final Multiome Object Creation)    |
+------------------------------------+
       |
       v
+------------------------------------+
|  BLOCK 6 & 7: Final QC & Recluster   |
|  (Interactive QC then Recluster)     |
+------------------------------------+
       |
       v
+------------------------------------+
|  BLOCK 8: Annotation & Subclustering | <-- Cell Typist, mapping_cell_types
|  (Genotype & Cell Identity Assign)   |
+------------------------------------+
       |
       v
+------------------------------------+
|  BLOCK 9: Genotype-Specific Export   | <-- Recluster on genotyped cells,
|  (Link peaks, Export h5ad for SIMBA)|     Export final files
+------------------------------------+
       |
       v
+------------------------------------+
|  BLOCK 10: Final RPCA & Graphs     | <-- 201_RPCA, Final_graphs.ipynb
+------------------------------------+
```

### scratch path /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed

########################## GEX libraries ###################################################################################################################
########################## GEX libraries ###################################################################################################################
########################## GEX libraries ###################################################################################################################

#### Index the barcode genome with cellranger

$ cellranger mkref --fasta /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa --genes /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/STAR.gtf --genome GFP_transgene_vCHEK2_and_DNMT3A_cellranger



#### 1. Modify targeted amplification of GEX to make it pass as CellRanger input

$ bash ~/Scripts/Wraper_scripts/127_cellranger_alignment_of_targeted_amp_GEX.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ targeted_amplicon_GEX /group/soranzo/manuel.tardaguila/Multiome/MCO_20250123/250124_A02059_0109_AHWTHYDSXC/adapter_trimmed_fastq/

$ bash ~/Scripts/Wraper_scripts/146_cellranger_alignment_of_targeted_amp_GEX_lymph.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ targeted_amplicon_GEX /group/soranzo/manuel.tardaguila/Multiome/QC_20250310/250310_A02059_0119_AH2N2TDMX2/adapter_trimmed_fastq/

$ bash ~/Scripts/Wraper_scripts/171_cellranger_alignment_of_targeted_amp_GEX_lymph.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ targeted_amplicon_GEX /group/soranzo/manuel.tardaguila/Multiome/RITM0029357/250414_A01481_0283_BHMKVYDRX5/fastq_raw/

#### 2. Run cellranger on the adapted reads

$ sbatch ~/Scripts/sbatch/18_Cell_ranger_count_for_CHEK2_targeted_amplicon_GEX_libraries_vMultiome.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01326 /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/cellranger/

$ sbatch ~/Scripts/sbatch/18_Cell_ranger_count_for_CHEK2_targeted_amplicon_GEX_libraries_vMultiome.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01327 /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/cellranger/

$ sbatch ~/Scripts/sbatch/18_Cell_ranger_count_for_CHEK2_targeted_amplicon_GEX_libraries_vMultiome.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01328 /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/cellranger/

$ sbatch ~/Scripts/sbatch/18_Cell_ranger_count_for_CHEK2_targeted_amplicon_GEX_libraries_vMultiome.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01329 /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/cellranger/

$ sbatch ~/Scripts/sbatch/18_Cell_ranger_count_for_CHEK2_targeted_amplicon_GEX_libraries_vMultiome.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01330 /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/cellranger/

$ sbatch ~/Scripts/sbatch/18_Cell_ranger_count_for_CHEK2_targeted_amplicon_GEX_libraries_vMultiome.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01331 /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/cellranger/

$ sbatch ~/Scripts/sbatch/18_Cell_ranger_count_for_CHEK2_targeted_amplicon_GEX_libraries_vMultiome.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01332 /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/cellranger/

$ sbatch ~/Scripts/sbatch/18_Cell_ranger_count_for_CHEK2_targeted_amplicon_GEX_libraries_vMultiome.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01333 /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/cellranger/



#### 3. Run CellBender to correct background of empty beads

# $ sbatch /home/manuel.tardaguila/Scripts/sbatch/15_cellbender_on_targeted_GEX_amplification.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01326

# $ sbatch /home/manuel.tardaguila/Scripts/sbatch/15_cellbender_on_targeted_GEX_amplification.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01327

# $ sbatch /home/manuel.tardaguila/Scripts/sbatch/15_cellbender_on_targeted_GEX_amplification.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01328

# $ sbatch /home/manuel.tardaguila/Scripts/sbatch/15_cellbender_on_targeted_GEX_amplification.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01329

# $ sbatch /home/manuel.tardaguila/Scripts/sbatch/15_cellbender_on_targeted_GEX_amplification.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01330

# $ sbatch /home/manuel.tardaguila/Scripts/sbatch/15_cellbender_on_targeted_GEX_amplification.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01331

# $ sbatch /home/manuel.tardaguila/Scripts/sbatch/15_cellbender_on_targeted_GEX_amplification.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01332

# $ sbatch /home/manuel.tardaguila/Scripts/sbatch/15_cellbender_on_targeted_GEX_amplification.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/targeted_amplicon_GEX/ MCO_01333


########################## GEX libraries ###################################################################################################################
########################## GEX libraries ###################################################################################################################
########################## GEX libraries ###################################################################################################################

### Block 1

$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01326
Submitted batch job 25021048
$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01327
Submitted batch job 25021049
$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01328
Submitted batch job 25021050
$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01329
Submitted batch job 25021051
$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01330
Submitted batch job 25021052
$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01331
Submitted batch job 25021053
$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01332
Submitted batch job 25021054
$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01333
Submitted batch job 25021055


### Block 2

$ bash ~/Scripts/Wraper_scripts/120_Seurat_first_pass_v3.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ processing_outputs

$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01326
$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01327
$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01328
$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01329
$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01330
$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01331
$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01332
$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ MCO_01333


$ bash ~/Scripts/Wraper_scripts/122_snATAC_pipeline.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ processing_outputs



$ sbatch ~/Scripts/sbatch/8_merge_atac_peaks_v3.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/


$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v4_MULTIOME.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa MCO_01326
$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v4_MULTIOME.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa MCO_01327
$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v4_MULTIOME.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa MCO_01328
$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v4_MULTIOME.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa MCO_01329
$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v4_MULTIOME.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa MCO_01330
$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v4_MULTIOME.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa MCO_01331
$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v4_MULTIOME.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa MCO_01332
$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v4_MULTIOME.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa MCO_01333

$ bash ~/Scripts/Wraper_scripts/119_Filter_Larry_and_graphs_v4.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/deconvolute_LARRY/ count_and_filter /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/deconvolute_LARRY/


### Block 3

$ bash ~/Scripts/Wraper_scripts/123_Amulet_v2.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ processing_outputs

### Block 4

$ bash ~/Scripts/Wraper_scripts/125_Seurat_second_pass_v3.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/

### Block 5

$ bash ~/Scripts/Wraper_scripts/126_merge_pre_merged_per_sample_v3.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ processing_outputs /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/processing_outputs/merged.atac_fragments.tsv.gz

### Block 6: Final QC

----> Jupyter notebook: Final_QC_in_the_merged_object.ipynb

### Block 7: recluster_after_QC

$ bash ~/Scripts/Wraper_scripts/200_Multiome_recluster_after_QC.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ processing_outputs

### Block 8 Cell_Typist, mapping_cell_types, Subclustering and barcode assignation

----> Jupyter notebook: Cell_Typist_triple_prediction_cell_identity.ipynb
----> Jupyter notebook: mapping_cell_types.ipynb
----> Jupyter notebook: Subclustering.ipynb
----> Jupyter notebook: notebook_to_assign_barcodes.ipynb

### Block 9 Recluster only on genotyped cells, call peaks using cell annotation and export final ATAC and RNA (corrected by Cell Bender but unormalized, the RNA layer) files to h5ad files for SIMBA

$ bash /home/manuel.tardaguila/Scripts/Wraper_scripts/202_Recluster_and_export_h5ad_after_G_v2.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ processing_outputs

### Block  10 Rpca and Final graphs

$ bash ~/Scripts/Wraper_scripts/201_RPCA_and_clustering_on_SCT.sh /scratch/manuel.tardaguila/2025_hESC_competition_assays_reanalysed/ processing_outputs

----> Jupyter notebook: Final_graphs.ipynb
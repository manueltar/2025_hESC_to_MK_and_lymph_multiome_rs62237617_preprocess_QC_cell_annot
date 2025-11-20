import os
import argparse
import sys
import pandas as pd
import anndata as ad
from scipy.sparse import csc_matrix, csr_matrix
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri 
from rpy2.rinterface_lib.embedded import RRuntimeError
import numpy as np 

# ----------------------------------------------------------------------

# --- Core Rpy2 Setup and Matrix Loading ---

def setup_r_environment():
    """Sets up the R environment necessary for rpy2 operations."""
    try:
        ro.r('library(Matrix)')
        ro.r('library(base)')
        print("R environment initialized with Matrix and base packages.")
    except RRuntimeError as e:
        print(f"FATAL R ERROR: Could not load R packages. Check R installation and libraries. Error: {e}", file=sys.stderr)
        sys.exit(1)

def load_r_matrix_to_csr(matrix_file: str, modality_name: str) -> tuple[csr_matrix, list]:
    """
    Loads an R dgCMatrix from an RDS file and converts it to a SciPy CSR matrix (Cells x Features).
    Returns the matrix and the feature names.
    """
    if not os.path.exists(matrix_file):
        raise FileNotFoundError(f"Matrix file for {modality_name} not found: {matrix_file}")

    with localconverter(ro.default_converter + pandas2ri.converter):
        print(f"Reading {modality_name} matrix from: {matrix_file}")
        try:
            r_matrix = ro.r(f'readRDS("{matrix_file}")')
        except RRuntimeError as e:
            raise RRuntimeError(f"R failed to read {matrix_file}. Error: {e}")

        R_inherits = ro.r['inherits']
        if not R_inherits(r_matrix, "dgCMatrix")[0]:
            raise TypeError(f"The loaded R {modality_name} matrix is not in dgCMatrix format.")

        # Extract components
        matrix_csr_data = r_matrix.do_slot("x")
        matrix_indices = r_matrix.do_slot("i")
        matrix_indptr = r_matrix.do_slot("p")
        matrix_shape = r_matrix.do_slot("Dim")

        # FIX 1: Ensure matrix data is float32 for H5AD compatibility
        csc_matrix_temp = csc_matrix(
            (matrix_csr_data.astype('float32'), matrix_indices, matrix_indptr),
            shape=tuple(matrix_shape)
        )
        X_matrix = csc_matrix_temp.T.tocsr()
        print(f"{modality_name} Matrix loaded successfully. Shape (Cells x Features): {X_matrix.shape}")

        R_rownames = ro.r['rownames']
        r_feature_names = list(R_rownames(r_matrix))
        
        # FIX 2 & 6: Robustly ensure all feature names are decoded Python 'str' objects
        feature_names = np.array([str(name).decode('utf-8') if isinstance(name, bytes) else str(name) 
                                  for name in r_feature_names], dtype='U').tolist()

    return X_matrix, feature_names

# ----------------------------------------------------------------------

# --- New Function: Peaks.bed Generation ---

def generate_peaks_bed(peak_features: list, output_dir: str):
    """
    Takes a list of peak strings (e.g., 'chr1-12345-67890') and writes a
    standard 5-column BED file (chrom, start, end, name, score).
    """
    print("Generating peaks.bed file...")

    # Split the peak strings (assuming format 'chr-start-end')
    peaks_df = pd.DataFrame({'peak_id': peak_features})
    peaks_df[['chrom', 'start', 'end']] = peaks_df['peak_id'].str.split('-', expand=True)

    # Convert start/end to integers
    peaks_df['start'] = pd.to_numeric(peaks_df['start'], errors='coerce').astype(int)
    
    # FIX 3: Convert 1-based start (from R) to 0-based (for BED)
    peaks_df['start'] = peaks_df['start'] - 1 

    peaks_df['end'] = pd.to_numeric(peaks_df['end'], errors='coerce').astype(int)

    # Create the final BED format: chrom, start, end, name, score (0)
    bed_output = peaks_df[['chrom', 'start', 'end']].copy()
    bed_output['name'] = peaks_df['peak_id']
    bed_output['score'] = 0

    bed_output_path = os.path.join(output_dir, "peaks.bed")

    # Write the file without header and using tab separator
    bed_output.to_csv(bed_output_path, sep='\t', header=False, index=False)
    print(f"Successfully exported peaks.bed to: {bed_output_path}")

# ----------------------------------------------------------------------

# --- Multi-Modal AnnData Assembly ---

def load_and_assemble_two_anndata(
    rna_matrix_file: str,
    atac_matrix_file: str,
    metadata_file: str
) -> tuple[ad.AnnData, ad.AnnData, list]:
    """
    Loads all data and returns two separate AnnData objects and the ATAC peak names.
    """
    # --- 1. Load RNA/SCT Data ---
    X_rna, feature_names_rna = load_r_matrix_to_csr(rna_matrix_file, "RNA/SCT")

    # --- 2. Load ATAC Data ---
    X_atac, feature_names_atac = load_r_matrix_to_csr(atac_matrix_file, "ATAC")

    # --- 3. Load Metadata (obs) ---
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

    with localconverter(ro.default_converter + pandas2ri.converter):
        print(f"Reading metadata from: {metadata_file}")
        r_metadata = ro.r(f'readRDS("{metadata_file}")')
        obs_df = ro.conversion.rpy2py(r_metadata)
        obs_df.index.name = 'cell_barcode'
        print(f"Metadata loaded successfully. Original Shape: {obs_df.shape}")
        
    # NEW FIX: Filter to ONLY include the desired columns
    required_columns = ['time_point', 'clone_line', 'Genotype', 'Construction_annotation', 'seurat_clusters']
    
    # Check for missing columns (in case R output changed)
    missing_cols = [col for col in required_columns if col not in obs_df.columns]
    if missing_cols:
        print(f"WARNING: The following required metadata columns were not found: {', '.join(missing_cols)}", file=sys.stderr)
        
    # Filter the DataFrame to only include the required columns that exist
    obs_df = obs_df[[col for col in required_columns if col in obs_df.columns]].copy()
    print(f"Metadata filtered. New Shape: {obs_df.shape} (Only requested columns retained)")


    # FIX 4 (Cleaned): Final Metadata Dtype Sanitization (Should be safer now)
    for col in list(obs_df.columns): 
        # Convert all factor/character (object) columns to the H5AD-compatible category type
        if obs_df[col].dtype == 'object' or isinstance(obs_df[col].dtype, pd.CategoricalDtype):
            obs_df[col] = obs_df[col].astype(str).astype('category')
            
    # FIX 7: Explicitly cast the shared cell barcode index to string
    obs_df.index = obs_df.index.astype(str)

    # --- 4. Validate Dimensions ---
    num_cells_rna = X_rna.shape[0]
    num_cells_atac = X_atac.shape[0]
    num_cells_meta = obs_df.shape[0]

    if not (num_cells_rna == num_cells_atac == num_cells_meta):
        raise ValueError(f"FATAL DIMENSION MISMATCH in cell counts: \n"
                         f"RNA cells: {num_cells_rna}, ATAC cells: {num_cells_atac}, Metadata cells: {num_cells_meta}")

    # --- 5. Create Feature Metadata ---
    # RNA: Use gene names as index 
    var_rna_df = pd.DataFrame(index=feature_names_rna)
    var_rna_df.index.name = 'gene_symbol'

    # FINAL WORKAROUND for ATAC serialization error:
    # Use a simple numeric index and store the complex peak regions as a column.
    var_atac_df = pd.DataFrame({
        'peak_region': np.array(feature_names_atac, dtype='U') 
    })
    # The index will be the default RangeIndex (0, 1, 2, ...), which is guaranteed to serialize.

    # --- 6. Final AnnData Assembly ---
    # Both AnnData objects get the *same* cleaned obs_df
    rna_adata = ad.AnnData(X=X_rna, obs=obs_df.copy(), var=var_rna_df)
    
    # ATAC: Use the var_atac_df with the simple numeric index
    atac_adata = ad.AnnData(X=X_atac, obs=obs_df.copy(), var=var_atac_df)

    return rna_adata, atac_adata, feature_names_atac

# ----------------------------------------------------------------------

def main():
    """Main function to parse arguments and execute the export."""
    parser = argparse.ArgumentParser(
        description="Converts R-exported sparse matrices (RNA/SCT and ATAC) and metadata RDS files into two separate H5AD files AND a peaks.bed file.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # --- INPUT/OUTPUT ARGUMENTS ---
    parser.add_argument(
        '--output-dir',
        type=str,
        default='/scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/processing_outputs/',
        help='Directory where the final H5AD files and peaks.bed will be saved.'
    )
    parser.add_argument(
        '--rna-matrix-file',
        type=str,
        required=True,
        help='FULL PATH to the RNA/SCT sparse matrix RDS file.'
    )
    parser.add_argument(
        '--atac-matrix-file',
        type=str,
        required=True,
        help='FULL PATH to the ATAC sparse matrix RDS file.'
    )
    parser.add_argument(
        '--metadata-file',
        type=str,
        required=True,
        help='FULL PATH to the cell metadata RDS file.'
    )
    parser.add_argument(
        '--base-output-name',
        type=str,
        default='merged_multiome',
        help='Base name for the output files (e.g., "base" will result in base_rna.h5ad and base_atac.h5ad).'
    )

    # --- HPC/RESOURCE ARGUMENTS (Kept) ---
    parser.add_argument('--cores', type=int, default=1, help='Number of cores/CPUs requested.')
    parser.add_argument('--memory', type=str, default='8G', help='Total memory allocated (e.g., "8G").')

    args = parser.parse_args()

    # --- Define Paths ---
    rna_matrix_path = args.rna_matrix_file
    atac_matrix_path = args.atac_matrix_file
    metadata_path = args.metadata_file

    rna_h5ad_output_path = os.path.join(args.output_dir, f"{args.base_output_name}_rna.h5ad")
    atac_h5ad_output_path = os.path.join(args.output_dir, f"{args.base_output_name}_atac.h5ad")
    output_dir = args.output_dir

    print(f"--- UNIFIED EXPORT STARTING ---")
    print(f"Output Directory: {output_dir}")
    print(f"RNA Output: {rna_h5ad_output_path}")
    print(f"ATAC Output: {atac_h5ad_output_path}")
    print("-" * 30)

    try:
        setup_r_environment()

        # Load and assemble two separate AnnData objects and get peaks
        rna_adata, atac_adata, feature_names_atac = load_and_assemble_two_anndata(
            rna_matrix_path,
            atac_matrix_path,
            metadata_path
        )

        # --- STEP 1: Generate Peaks.bed file ---
        generate_peaks_bed(feature_names_atac, output_dir)

        # --- STEP 2: Save RNA AnnData ---
        print(f"Saving RNA AnnData object to: {rna_h5ad_output_path}")
        rna_adata.write_h5ad(rna_h5ad_output_path)

        # --- STEP 3: Save ATAC AnnData ---
        print(f"Saving ATAC AnnData object to: {atac_h5ad_output_path}")
        atac_adata.write_h5ad(atac_h5ad_output_path)

        print("\n---------------------------------------------------------")
        print("EXPORT SUCCESS! 🎉 All files generated.")
        print("---------------------------------------------------------")

    except Exception as e:
        print(f"\nCRITICAL SCRIPT FAILURE: {e}", file=sys.stderr)
        print("Export failed. Please check file paths, R/Python libraries, and matrix dimensions.", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()

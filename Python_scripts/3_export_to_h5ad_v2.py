import os
import argparse
import sys
import pandas as pd
import anndata as ad
from scipy.sparse import csc_matrix
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri
from rpy2.rinterface_lib.embedded import RRuntimeError

def setup_r_environment():
    """Sets up the R environment necessary for rpy2 operations."""
    try:
        # Load required R packages
        ro.r('library(Matrix)')
        ro.r('library(base)')
        print("R environment initialized with Matrix and base packages.")
    except RRuntimeError as e:
        print(f"FATAL R ERROR: Could not load R packages. Check R installation and libraries. Error: {e}", file=sys.stderr)
        sys.exit(1)

def load_and_assemble_data(matrix_file: str, metadata_file: str) -> ad.AnnData:
    """
    Loads R sparse matrix and metadata, converts them to AnnData-compatible formats,
    and assembles the final AnnData object.
    """
    if not all(os.path.exists(f) for f in [matrix_file, metadata_file]):
        raise FileNotFoundError("One or both input RDS files were not found.")

    # Use the modern rpy2 conversion context manager
    with localconverter(ro.default_converter + pandas2ri.converter):
        print(f"Reading matrix from: {matrix_file}")
        r_matrix = ro.r(f'readRDS("{matrix_file}")')

        # --- 1. Validate and Extract Matrix Components ---
        R_inherits = ro.r['inherits']
        if not R_inherits(r_matrix, "dgCMatrix")[0]:
            raise TypeError("The loaded R matrix is not in dgCMatrix format.")

        matrix_csr = r_matrix.do_slot("x")
        matrix_indices = r_matrix.do_slot("i")
        matrix_indptr = r_matrix.do_slot("p")
        matrix_shape = r_matrix.do_slot("Dim")

        # CRITICAL FIX: R's dgCMatrix is CSC. Construct as CSC, then transpose to CSR (NxG).
        csc_X_matrix = csc_matrix(
            (matrix_csr, matrix_indices, matrix_indptr),
            shape=tuple(matrix_shape)
        )
        X_matrix = csc_X_matrix.T
        print(f"Matrix loaded successfully. Shape (Cells x Genes): {X_matrix.shape}")

        # --- 2. Load Metadata (obs) ---
        print(f"Reading metadata from: {metadata_file}")
        r_metadata = ro.r(f'readRDS("{metadata_file}")')
        obs_df = ro.conversion.rpy2py(r_metadata)
        obs_df.index.name = 'cell_barcode'
        print(f"Metadata loaded successfully. Shape: {obs_df.shape}")

        # --- 3. Create Gene Metadata (var) ---
        R_rownames = ro.r['rownames']
        r_feature_names = R_rownames(r_matrix)
        var_df = pd.DataFrame(index=list(r_feature_names))
        print(f"Feature metadata created. Shape: {var_df.shape}")

    # --- 4. Final AnnData Assembly ---
    if X_matrix.shape[0] != obs_df.shape[0]:
        raise ValueError("FATAL DIMENSION MISMATCH: Number of cells in matrix does not match metadata rows.")
    if X_matrix.shape[1] != var_df.shape[0]:
        raise ValueError("FATAL DIMENSION MISMATCH: Number of genes in matrix does not match feature rows.")

    adata = ad.AnnData(
        X=X_matrix,
        obs=obs_df,
        var=var_df
    )
    return adata

def main():
    """Main function to parse arguments and execute the export."""
    # Use argparse to make the script flexible and runnable
    parser = argparse.ArgumentParser(
        description="Converts R-exported sparse matrix and metadata RDS files into a single H5AD file.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # --- INPUT/OUTPUT ARGUMENTS ---
    parser.add_argument(
        '--output-dir',
        type=str,
        default='/scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/processing_outputs/',
        help='Directory where the final H5AD file will be saved.'
    )
    # MODIFICATION 1: Input files now require full paths
    parser.add_argument(
        '--matrix-file',
        type=str,
        required=True, # Made required since it's an absolute path
        help='FULL PATH to the sparse matrix RDS file.'
    )
    # MODIFICATION 1: Input files now require full paths
    parser.add_argument(
        '--metadata-file',
        type=str,
        required=True, # Made required since it's an absolute path
        help='FULL PATH to the cell metadata RDS file.'
    )
    parser.add_argument(
        '--output-name',
        type=str,
        default='merged_sct_final_export.h5ad',
        help='Name for the final H5AD output file.'
    )
    
    # --- HPC/RESOURCE ARGUMENTS (New) ---
    parser.add_argument(
        '--cores',
        type=int,
        default=1,
        help='Number of cores/CPUs requested for this job (for documentation/potential downstream use).'
    )
    parser.add_argument(
        '--memory',
        type=str,
        default='8G',
        help='Total memory allocated for this job (e.g., "8G", "16GB").'
    )


    args = parser.parse_args()

    # --- Define Paths (Matrix/Metadata are already full paths) ---
    matrix_path = args.matrix_file
    metadata_path = args.metadata_file
    # Output path is constructed from the output directory and the output name
    h5ad_output_path = os.path.join(args.output_dir, args.output_name)

    print(f"--- H5AD EXPORT STARTING ---")
    print(f"Output Directory: {args.output_dir}")
    print(f"Input Matrix File: {matrix_path}")
    print(f"Input Metadata File: {metadata_path}")
    print(f"Resource Request: {args.cores} cores, {args.memory} memory")
    print("-" * 30)

    try:
        setup_r_environment()

        adata = load_and_assemble_data(matrix_path, metadata_path)

        # Save to H5AD
        print(f"Saving final AnnData object to: {h5ad_output_path}")
        adata.write_h5ad(h5ad_output_path)

        print("\n---------------------------------------------------------")
        print("EXPORT SUCCESS! 🎉")
        print(f"File: {h5ad_output_path}")
        print(f"AnnData Summary: {adata}")
        print("---------------------------------------------------------")

    except Exception as e:
        print(f"\nCRITICAL SCRIPT FAILURE: {e}", file=sys.stderr)
        print("Export failed. Please check file paths and library installations.", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()

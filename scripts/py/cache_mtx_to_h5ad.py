# Save an anndata object from market matrix input. Also adds spatial coordinates to the object

import pandas as pd
import argparse
from scanpy import read_mtx, pl
from numpy import intersect1d
import matplotlib.pyplot as plt

# Example usage:
""" 
python scripts/py/cache_mtx_to_h5ad.py \
    --mat_in input_matrix.mtx \
    --feat_in input_features.tsv \
    --bc_in input_barcodes.txt \
    --bc_map input_spatial_map.tsv \
    --ad_out output_anndata.h5ad \
    --feat_col 1 \
    --remove_zero_features
"""

def plot_qc_metrics(adata, output_file):
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot total counts per cell (knee plot)
    sorted_counts = adata.obs['n_counts'].sort_values(ascending=False)
    axes[0, 0].plot(range(len(sorted_counts)), sorted_counts)
    axes[0, 0].set_xlabel('Cell Rank')
    axes[0, 0].set_ylabel('Total Counts')
    axes[0, 0].set_title('Knee Plot of Total Counts per Cell')

    # Plot number of genes per cell (knee plot)
    sorted_genes = adata.obs['n_genes'].sort_values(ascending=False)
    axes[0, 1].plot(range(len(sorted_genes)), sorted_genes)
    axes[0, 1].set_xlabel('Cell Rank')
    axes[0, 1].set_ylabel('Number of Genes')
    axes[0, 1].set_title('Knee Plot of Number of Genes per Cell')

    # Create spatial plots using scanpy
    pl.embedding(adata, basis='spatial', color='n_counts', ax=axes[1, 0], show=False)
    axes[1, 0].set_title('Spatial Map of Total Counts')

    pl.embedding(adata, basis='spatial', color='n_genes', ax=axes[1, 1], show=False)
    axes[1, 1].set_title('Spatial Map of Number of Genes')

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"QC plots saved to {output_file}")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Process spatial transcriptomics data."
    )
    parser.add_argument(
        "--mat_in", required=True, help="Input count matrix file (mtx format)"
    )
    parser.add_argument(
        "--feat_in", required=True, help="Input feature file (tsv format)"
    )
    parser.add_argument(
        "--bc_in", required=True, help="Input barcode file (txt format)"
    )
    parser.add_argument(
        "--bc_map", required=True, help="Input spatial map file (tsv format)"
    )
    parser.add_argument(
        "--ad_out", required=True, help="Output AnnData file (h5ad format)"
    )
    parser.add_argument(
        "--feat_cols",
        type=int,
        nargs="+",
        default=[1],
        help="Feature column indices in the feature file (default: [1])",
    )
    parser.add_argument(
        "--transpose",
        type=bool,
        default=False,
        help="Transpose count matrix? (default: False)",
    )
    parser.add_argument(
        "--remove_zero_features",
        action="store_true",
        help="Remove observations with zero features detected (default: False)",
    )
    parser.add_argument(
        "--plot_qc",
        action="store_true",
        help="Plot QC metrics and save the plots (default: False)",
    )
    parser.add_argument(
        "--qc_plot_file",
        type=str,
        default=None,
        help="Filename for the QC plots (default: {output_dir}/qc_plots.png)",
    )
    return parser.parse_args()

def main(
    mat_in,
    feat_in,
    bc_in,
    bc_map,
    ad_out,
    feat_cols=[1],
    transpose=True,
    remove_zero_features=False,
    verbose=True,
    plot_qc=False,
    qc_plot_file=None,
):
    # Count matrix
    adata = read_mtx(mat_in)

    # Transpose for STAR inputs...
    if transpose:
        print("transposing count matrix...")
        print("")
        adata = adata.transpose()

    # Features
    feature_df = pd.read_csv(feat_in, sep="\t", header=None)

    # Set the first feature column as var_names
    adata.var_names = feature_df.iloc[:, feat_cols[0]].values

    # Add all specified feature columns to adata.var
    for i, col in enumerate(feat_cols):
        col_name = f"feature_col_{col}"
        adata.var[col_name] = feature_df.iloc[:, col].values

    # Barcodes
    adata.obs_names = pd.read_csv(bc_in, sep="\t", header=None)[0].tolist()

    # Add spatial location
    spatial_data = pd.read_csv(
        bc_map, sep="\t", header=None, names=["barcode", "x", "y"]
    )

    if verbose:
        print(
            f"Found {len(adata.var_names)} features in count matrix.\n"
            f"Found {len(adata.obs_names)} cell barcodes in count matrix.\n"
            f"Found {spatial_data.shape[0]} cell barcodes in whitelist/map.\n"
        )

    # Set the cell barcode as index
    spatial_data.set_index("barcode", inplace=True)

    # Check if the barcodes in the AnnData object match the ones in the spatial data
    if not all(adata.obs_names.isin(spatial_data.index)):
        print(
            "Warning: Not all barcodes in the AnnData object match the ones in the barcode map..."
        )

    # Add the spatial coordinates to the AnnData object
    common_labels = intersect1d(adata.obs_names, spatial_data.index)
    adata = adata[common_labels, :]
    spatial_coord = spatial_data.loc[common_labels, ["x", "y"]]

    print(
        f"{len(common_labels)} / {len(adata.obs_names)} in output AnnData object found in barcode map..."
    )

    adata.obsm["spatial"] = spatial_coord.to_numpy()

    # Remove observations with zero features detected
    if remove_zero_features:
        adata = adata[adata.X.sum(axis=1) > 0, :]

    # Calculate QC metrics
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)

    # Plot QC metrics if requested
    if plot_qc:
        if qc_plot_file is None:
            output_dir = "/".join(ad_out.split("/")[:-1])
            qc_plot_file = f"{output_dir}/qc_plots.png"
        plot_qc_metrics(adata, qc_plot_file)

    # Write output
    print(f"Writing to {ad_out}")
    adata.write(ad_out)


if __name__ == "__main__":
    args = parse_args()
    print(
        f"Matrix file:                  {args.mat_in}\n"
        f"Features/genes file:          {args.feat_in}\n"
        f"Barcodes file:                {args.bc_in}\n"
        f"Barcode map file:             {args.bc_map}\n"
        f"Output AnnData file:          {args.ad_out}\n"
        f"Feature column indices:       {args.feat_cols}\n"
        f"Transpose matrix:             {args.transpose}\n"
        f"Remove undetected features:   {args.remove_zero_features}\n"
        f"Plot QC metrics:              {args.plot_qc}\n"
        f"QC plot file:                 {args.qc_plot_file}\n"
    )
    main(
        args.mat_in,
        args.feat_in,
        args.bc_in,
        args.bc_map,
        args.ad_out,
        args.feat_cols,
        args.transpose,
        args.remove_zero_features,
        plot_qc=args.plot_qc,
        qc_plot_file=args.qc_plot_file,
    )

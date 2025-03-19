# Save an anndata object from market matrix input. Also adds spatial coordinates to the object

import pandas as pd
import argparse
import gzip
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

    # default = 120000 / num_points
    num_points = adata.shape[0]
    if num_points > 5000:
        point_size = 150000 / num_points
    else:
        point_size = max(1, 150000 / num_points)

    # Plot total counts per cell (knee plot)
    sorted_counts = adata.obs["n_counts"].sort_values(ascending=False)
    axes[0, 0].plot(range(len(sorted_counts)), sorted_counts)
    axes[0, 0].set_xlabel("Cell Rank")
    axes[0, 0].set_ylabel("Total Counts")
    axes[0, 0].set_title("Knee Plot of Total Counts per Cell")

    # Plot number of genes per cell (knee plot)
    sorted_genes = adata.obs["n_genes"].sort_values(ascending=False)
    axes[0, 1].plot(range(len(sorted_genes)), sorted_genes)
    axes[0, 1].set_xlabel("Cell Rank")
    axes[0, 1].set_ylabel("Number of Genes")
    axes[0, 1].set_title("Knee Plot of Number of Genes per Cell")

    # Create spatial plots using scanpy
    pl.embedding(
        adata,
        basis="spatial",
        color="n_counts",
        ax=axes[1, 0],
        show=False,
        size=point_size,
    )
    axes[1, 0].set_title("Spatial Map of Total Counts")

    pl.embedding(
        adata,
        basis="spatial",
        color="n_genes",
        ax=axes[1, 1],
        show=False,
        size=point_size,
    )
    axes[1, 1].set_title("Spatial Map of Number of Genes")

    plt.tight_layout()
    plt.savefig(output_file, dpi=400)
    plt.close()
    print(f"QC plots saved to {output_file}")


def load_gtf_to_dataframe(gtf_file, feature_type="all", seqname_filter="all", unwrap_attributes=True):
    """
    Load a GTF file (including gzipped GTF files) into a Pandas DataFrame.

    Parameters:
        gtf_file (str): Path to the GTF file. Can be a plain text or gzipped file.
        feature_type (str): Type of feature to load (e.g., "gene", "exon"). 
                            Default is "all", which loads all features.
        seqname_filter (str): Sequence name to filter by (e.g., "chr1", "chr2"). 
                              Default is "all", which loads all sequence names.
        unwrap_attributes (bool): If True, unpacks the attributes column into separate columns.
                                  Default is True.

    Returns:
        pd.DataFrame: A DataFrame containing the GTF data.
    """
    def parse_attributes(attributes_str):
        """
        Parse the attributes column of a GTF file into a dictionary.
        """
        attributes = {}
        for attribute in attributes_str.strip().split(';'):
            if attribute.strip():
                key, value = attribute.strip().split(' ', 1)
                attributes[key] = value.strip('"')
        return attributes

    # Open the file (support for gzipped files)
    open_func = gzip.open if gtf_file.endswith('.gz') else open
    with open_func(gtf_file, 'rt') as f:
        # Read the GTF file into a DataFrame
        gtf_data = []
        attributes_data = []
        for line in f:
            if line.startswith('#'):
                continue  # Skip comment lines
            fields = line.strip().split('\t')
            if feature_type != "all" and fields[2] != feature_type:
                continue  # Skip features that don't match the specified type
            if seqname_filter != "all" and fields[0] != seqname_filter:
                continue  # Skip sequence names that don't match the filter
            attributes = parse_attributes(fields[8])
            gtf_data.append(fields[:8])
            attributes_data.append(attributes)
        
        # Create the main DataFrame
        columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']
        gtf_df = pd.DataFrame(gtf_data, columns=columns)
        
        if unwrap_attributes:
            # Create a DataFrame for attributes and concatenate with the main DataFrame
            attributes_df = pd.DataFrame(attributes_data)
            gtf_df = pd.concat([gtf_df, attributes_df], axis=1)
    
    return gtf_df

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
    parser.add_argument(
        "--gtf_file",
        type=str,
        default=None,
        help="Path to GTF file (optional)."
    )
    parser.add_argument(
        "--gtf_feature_type",
        type=str,
        default="gene",
        help="Feature type to load from GTF (default: gene)."
    )
    parser.add_argument(
        "--gtf_id",
        type=str,
        default="gene_name",
        help="Column in GTF used to match var_names"
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
    gtf_file=None,
    gtf_feature_type="gene",
    gtf_id="gene_name"
):
    # Count matrix
    adata = read_mtx(mat_in)

    # Transpose for STAR inputs...
    if transpose:
        # Convert to CSC format to avoid memory errors when transposing
        print("Converting count matrix to CSC & transposing...")
        adata.X = adata.X.tocsc()
        adata = adata.transpose()
        print("    Done.")

    # Save raw counts
    adata.raw = adata

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
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1)

    # Plot QC metrics if requested
    if plot_qc:
        if qc_plot_file is None:
            output_dir = "/".join(ad_out.split("/")[:-1])
            qc_plot_file = f"{output_dir}/qc_plots.png"
        plot_qc_metrics(adata, qc_plot_file)

    # Add feature info from gtf if provided
    if gtf_file:
        print(f"Loading GTF annotations from {gtf_file}...")
        gtf_df = load_gtf_to_dataframe(
            gtf_file,
            feature_type=gtf_feature_type
        )
        if gtf_id in gtf_df.columns:
            print(f"Reordering GTF rows using '{gtf_id}'...")
            gtf_df_ordered = gtf_df.set_index(gtf_id).reindex(adata.var_names)
            for col in gtf_df_ordered.columns:
                adata.var[col] = gtf_df_ordered[col].astype(str)
        else:
            print(f"Warning: '{gtf_id}' not found in GTF columns.")

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
        f"GTF file:                     {args.gtf_file}\n"
        f"GTF feature type:             {args.gtf_feature_type}\n"
        f"GTF ID column:                {args.gtf_id}\n"
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
        gtf_file=args.gtf_file,
        gtf_feature_type=args.gtf_feature_type,
        gtf_id=args.gtf_id
    )

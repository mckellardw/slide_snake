# Save an anndata object from market matrix input. Also adds spatial coordinates to the object

# import sys
# import gzip
import pandas as pd
import argparse
from scanpy import read_mtx
from numpy import intersect1d

# Example usage:
""" 
python script.py \
    --mat_in input_matrix.mtx \
    --feat_in input_features.tsv \
    --bc_in input_barcodes.txt \
    --bc_map input_spatial_map.tsv \
    --ad_out output_anndata.h5ad \
    --feat_col 1 \
    --remove_zero_features
"""


def main(
    mat_in,
    feat_in,
    bc_in,
    bc_map,
    ad_out,
    feat_col=1,
    remove_zero_features=False,
    verbose=True,
):
    # Count matrix
    adata = read_mtx(mat_in)

    # Transpose for STAR inputs...
    if "Solo.out" in mat_in:
        print("transposing count matrix...")
        adata = adata.transpose()

    # Features
    # adata.var_names = pd.read_csv(feat_in, sep="\t", header=None, usecols=[feat_col], squeeze=True).values    # pandas v1.#.#
    adata.var_names = (
        pd.read_csv(feat_in, sep="\t", header=None, usecols=[feat_col]).squeeze().values
    )  # pandas v2.#.#

    # Barcodes
    # adata.obs_names = pd.read_csv(bc_in, sep="\t", header=None, squeeze=True).values    # pandas v1.#.#
    adata.obs_names = pd.read_csv(bc_in, sep="\t", header=None)[
        0
    ].tolist()  # pandas v2.#.#

    # Add spatial location
    spatial_data = pd.read_csv(
        bc_map, sep="\t", header=None, names=["barcode", "x", "y"]
    )

    if verbose:
        print(f"Found {spatial_data.shape[0]} cell barcodes in whitelist/map.")

    # Set the cell barcode as index
    spatial_data.set_index("barcode", inplace=True)

    # Check if the barcodes in the AnnData object match the ones in the spatial data
    if not all(adata.obs_names.isin(spatial_data.index)):
        print(
            "Warning: Not all barcodes in the AnnData object match the ones in the spatial data."
        )

    # Add the spatial coordinates to the AnnData object
    common_labels = intersect1d(adata.obs_names, spatial_data.index)
    adata = adata[common_labels, :]
    spatial_coord = spatial_data.loc[common_labels, ["x", "y"]]

    print(f"{len(common_labels)} / {len(adata.obs_names)} found in barcode map...")

    # spatial_coord = spatial_data.loc[adata.obs_names, ['x', 'y']] #w/ out filtering for intersection
    adata.obsm["spatial"] = spatial_coord.to_numpy()

    # Remove observations with zero features detected
    if remove_zero_features:
        adata = adata[adata.X.sum(axis=1) > 0, :]

    # Write output
    adata.write(ad_out)


if __name__ == "__main__":
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
        "--feat_col",
        type=int,
        default=1,
        help="Feature column index in the feature file (default: 1)",
    )
    parser.add_argument(
        "--remove_zero_features",
        action="store_true",
        help="Remove observations with zero features detected (default: False)",
    )

    args = parser.parse_args()
    print(
        f"Matrix file:                  {args.mat_in}\n"
        f"Features/genes file:          {args.feat_in}\n"
        f"Barcodes file:                {args.bc_in}\n"
        f"Barcode map file:             {args.bc_map}\n"
        f"Output AnnData file:          {args.ad_out}\n"
        f"Feature column index:         {args.feat_col}\n"
        f"Remove undetected features:   {args.remove_zero_features}\n"
    )
    main(
        args.mat_in,
        args.feat_in,
        args.bc_in,
        args.bc_map,
        args.ad_out,
        args.feat_col,
        args.remove_zero_features,
    )

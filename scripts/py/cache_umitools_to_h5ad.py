# Save an anndata object from umi_tools count matrix input. Also adds spatial coordinates to the object

import pandas as pd
import argparse
from anndata import read_umi_tools
from numpy import intersect1d

# Example usage:
"""
python scripts/py/cache_umitools_to_h5ad.py \
    --mat_in umitools_counts.tsv \
    --bc_map input_spatial_map.tsv \
    --ad_out output_anndata.h5ad \
    --remove_zero_features
"""


def parse_args():
    parser = argparse.ArgumentParser(
        description="Process spatial transcriptomics data from umi_tools count matrix."
    )
    parser.add_argument(
        "--mat_in",
        required=True,
        help="Input count matrix file (tsv format from umi_tools)",
    )
    parser.add_argument(
        "--bc_map", required=True, help="Input spatial map file (tsv format)"
    )
    parser.add_argument(
        "--ad_out", required=True, help="Output AnnData file (h5ad format)"
    )
    parser.add_argument(
        "--remove_zero_features",
        action="store_true",
        help="Remove observations with zero features detected (default: False)",
    )
    return parser.parse_args()


def main(
    mat_in,
    bc_map,
    ad_out,
    remove_zero_features=False,
    verbose=True,
):
    # Count matrix
    adata = read_umi_tools(mat_in)

    # Add spatial location
    spatial_data = pd.read_csv(
        bc_map, sep="\t", header=None, names=["barcode", "x", "y"]
    )

    if verbose:
        print(
            f"Barcodes in count matrix:     {len(adata.obs_names):,}\n"
            f"Barcodes in spatial map:      {spatial_data.shape[0]:,}\n"
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

    # Write output
    print(f"Writing to {ad_out}")
    adata.write(ad_out)


if __name__ == "__main__":
    args = parse_args()

    verbose = True

    print(
        f"Matrix file:                  {args.mat_in}\n"
        f"Barcode map file:             {args.bc_map}\n"
        f"Output AnnData file:          {args.ad_out}\n"
        f"Remove undetected features:   {args.remove_zero_features}\n"
    )
    main(
        mat_in=args.mat_in,
        bc_map=args.bc_map,
        ad_out=args.ad_out,
        remove_zero_features=args.remove_zero_features,
        verbose=verbose,
    )

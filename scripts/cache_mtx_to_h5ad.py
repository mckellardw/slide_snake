import sys
from scanpy import read_mtx
import gzip
import pandas as pd

args = sys.argv
mat_in = args[1]
feat_in = args[2]
bc_in = args[3]
bb_map = args[4]
ad_out = args[5]
# var_names = args[3]

# adata = read_10x_mtx(
#     path=mat_in,
#     var_names=var_names,
#     make_unique=True,
#     cache=False
# )

# Count matrix
if mat_in.endswith('.gz'):
    adata = read_mtx(mat_in)
else:
    adata = read_mtx(gzip.open(mat_in, "rt"))

# Features
if feat_in.endswith('.gz'):
    with gzip.open(feat_in, "rt") as f:
        adata.var_names = [line.strip() for line in f]
else:
    adata.var_names = [line.strip() for line in open(feat_in)]

# Barcodes
if bc_in.endswith('.gz'):
    with gzip.open(bc_in, "rt") as f:
        adata.obs_names = [line.strip() for line in f]
else:
    adata.obs_names = [line.strip() for line in open(bc_in)]

# Add spatial location
spatial_data = pd.read_csv(
    bb_map, 
    sep="\t", 
    header=None, 
    names=["barcode", "x", "y"]
)

# Set the cell barcode as index
spatial_data.set_index("barcode", inplace=True)

# Check if the barcodes in the AnnData object match the ones in the spatial data
if not all(adata.obs_names.isin(spatial_data.index)):
    print("Warning: Not all barcodes in the AnnData object match the ones in the spatial data.")

# Add the spatial coordinates to the AnnData object
spatial_coord = spatial_data.loc[adata.obs_names, ['x', 'y']]
adata.obsm['spatial'] = spatial_coord.to_numpy()

# Write output
adata.write(ad_out)

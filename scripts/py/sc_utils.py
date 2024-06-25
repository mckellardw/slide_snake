# Utility functions for use with scanpy
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib.path import Path
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from typing import Union
import scipy.sparse as sp
from scipy.io import mmread

import scipy.sparse as sp
import gzip


#  Get # PCs for a given % of variance
def npcs(ADATA, var_perc=0.95, reduction="pca"):
    """
    Calculate the number of Principal Components (PCs) that contain a certain proportion of the variance.

    Parameters:
    -------
    ADATA -- An AnnData object. Must have already been processed with PCA and contain a 'pca' entry in its 'obsm' field.
    var_perc -- A float indicating the proportion of variance to be covered by the selected PCs. Default is 0.95.
    reduction -- A string indicating the type of dimensionality reduction to use. Default is 'pca'.

    Returns:
    -------
    n_pcs -- The number of PCs needed to cover the specified proportion of the variance. If the specified 'reduction' is not found, returns None.
    """
    from numpy import sum, var

    get_var = lambda i: var(ADATA.obsm[reduction][:, i])

    if ADATA.obsm[reduction] is None:
        print(f"Reduction '{reduction}', not found!")
        return None
    else:
        var_tmp = [get_var(i) for i in list(range(0, ADATA.obsm[reduction].shape[1]))]
        var_cut = var_perc * sum(var_tmp)
        n_pcs = 0
        var_sum = 0
        while var_sum < var_cut and n_pcs < ADATA.obsm[reduction].shape[1] - 1:
            var_sum = var_sum + var_tmp[n_pcs]
            n_pcs = n_pcs + 1

        return n_pcs


# Pattern matching for gene names.
## Returns the genes sorted by expression (high to low) by default
def grep_genes(
    adata,
    pattern="",
    assay=None,
    filter_pattern=None,
    sort_by="expression",
    verbose=True,
):
    """
    Search for genes in an AnnData object that match a given pattern.

    Parameters:
    adata (anndata.AnnData): The AnnData object containing the gene expression data.
    pattern (str or list of str, optional): The pattern or patterns to search for in the gene names. Default is "".
    assay (str, optional): The assay to use. Default is None, which uses "counts".
    filter_pattern (str or list of str, optional): A pattern or patterns to filter out from the gene names. Default is None.
    sort_by (str, optional): How to sort the output genes. Can be "expression" (sort by total expression across all cells) or "abc" (sort alphabetically). Default is "expression".
    verbose (bool, optional): Whether to print progress messages. Default is True.

    Returns:
    list of str: The names of the genes that match the pattern and do not match the filter_pattern, sorted according to sort_by.
    """
    if pattern == "":
        if verbose:
            print("Need a pattern to look for!")
        return None

    if assay is None:
        assay = "counts"

    genes = adata.var_names

    if verbose:
        print(f"Found {len(genes)} features in the assay '{assay}'...")

    if isinstance(pattern, list):
        if verbose:
            print("Looking for multiple patterns in these features...")
        out_genes = [gene for gene in genes if any(pat in gene for pat in pattern)]
    else:
        if verbose:
            print(f"Looking for '{pattern}' in these features...")
        out_genes = [gene for gene in genes if pattern in gene]

    if len(out_genes) == 0:
        print("Nothing found!\n")
        return None

    if filter_pattern is not None:
        if verbose:
            print(f"Removing features containing '{filter_pattern}' from output...\n")
        out_genes = [
            gene for gene in out_genes if not any(pat in gene for pat in filter_pattern)
        ]

    if sort_by == "expression":
        if verbose:
            print("Output genes sorted by expression...")
        gene_sums = np.asarray(adata.X.sum(axis=0)).squeeze()
        out_genes = sorted(
            out_genes, key=lambda gene: -gene_sums[adata.var_names == gene]
        )
    elif sort_by == "abc":
        if verbose:
            print("Output genes sorted alphabetically...")
        out_genes = sorted(out_genes)

    return out_genes


# Reorder a reduction by decreasing % variance
def reorder_reduction(ADATA, reduction="pca", verbose=False):
    """
    Re-order a dimensions of a reduction by decreasing % variance.

    Parameters:
    -------
    ADATA -- An AnnData object. Must have already been processed with PCA and contain a 'pca' entry in its 'obsm' field.
    reduction -- A string indicating the type of dimensionality reduction to use. Default is 'pca'.
    verbose -- A boolean to indicate whether to print the variance for each dimension. Default is False.

    This function doesn't return anything, but it modifies the AnnData object in place, re-ordering the dimensions
    of the specified reduction in the 'obsm' field based on their variance (in decreasing order).
    """
    from numpy import var, argsort

    if reduction in ADATA.obsm:
        get_var = lambda i: var(ADATA.obsm[reduction][:, i])
        var_tmp = [get_var(i) for i in list(range(0, ADATA.obsm[reduction].shape[1]))]
        if verbose:
            print("Reduction variance by dimension:")
            print(var_tmp)

        pc_order = argsort(var_tmp)[::-1]
        ADATA.obsm[reduction] = ADATA.obsm[reduction][:, pc_order]
    else:
        print(f"The reduction '{reduction}' was not found...")


# Read in a list of gene lists from .csv (each column is a gene list)
def read_csv_to_dict(filename, names2check=""):
    """
    Read in a list of gene lists from .csv (each column is a gene list).

    Parameters:
    -------
    filename -- A string specifying the location of the csv file.
    names2check -- A list of gene names to filter the dictionary by. If not specified, all gene names are included.

    Returns:
    -------
    dict_out -- A dictionary with column headers from the csv file as keys, and lists of genes as values. If names2check is specified, only genes in names2check are included in the lists.
    """
    import csv

    # Open the CSV file
    with open(filename, "r") as file:
        # Create a CSV reader object
        reader = csv.reader(file)

        # Read the first row as header
        header = next(reader)

        # Create an empty dictionary to store the columns
        dict_out = {col: [] for col in header}

        # Loop through each row in the CSV file
        for row in reader:
            # Loop through each column in the row
            for col, value in zip(header, row):
                # Add the value to the corresponding column in the dictionary
                if value:  # skip empty strings
                    dict_out[col].append(value)

    # Filter out unwanted entries based on the list in `names2check`
    if len(names2check) > 1:
        for KEY in dict_out.keys():
            dict_out[KEY] = [k for k in dict_out[KEY] if k in names2check]

    # Return the dictionary
    return dict_out


# Function to export DGEA results to a .csv file
def export_dgea_to_csv(
    adata: ad.AnnData, dgea_name, n_features, csv_out, axis=0, wide=False
):
    """
    Function to export DGEA (Differential Gene Expression Analysis) results to a .csv file.

    Parameters:
    -------
    adata -- An AnnData object which stores the gene expression data and metadata.
    dgea_name -- A string specifying the name of the DGEA results in adata.uns[].
    n_features -- An integer specifying the number of top features to be included in the exported .csv file.
    csv_out -- A string specifying the path to the output .csv file.
    axis -- An integer specifying how to write results for each group. If 1, results are written horizontally. If 0, results are written vertically. Default is 0.
    wide -- A boolean specifying the format of the .csv file. If True, the .csv file will have one row per feature, and each column will be a group. If False, the .csv file will have one row per group-feature combination. Default is False.

    Returns:
    -------
    The function doesn't return anything but writes the DGEA results to a .csv file.
    """
    import pandas as pd
    import scanpy as sc

    result = adata.uns[dgea_name]
    groups = result["names"].dtype.names

    if wide:
        celltype_markers = pd.DataFrame(
            {
                group + "_" + key[:-1]: result[key][group]
                for group in groups
                for key in ["names", "logfoldchanges", "pvals"]
            }
        ).head(n_features)
        celltype_markers.to_csv(csv_out, index=False)
    else:
        marker_list = list()
        for group in adata.uns[dgea_name]["names"].dtype.names:
            markers = sc.get.rank_genes_groups_df(
                adata, key=dgea_name, group=group
            ).head(n_features)
            markers["celltypes"] = group
            marker_list.append(markers)

        celltype_markers = pd.concat(marker_list, axis=axis)
        celltype_markers.to_csv(csv_out, index=False)


# Function to convert feature names
def convert_feature_names(
    adata: ad.AnnData,
    gtf_info: pd.DataFrame,
    from_col: str = "GENEID",
    to_col: str = "GeneSymbol",
    inplace: bool = True,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Function to convert feature names in an AnnData object using mapping provided in a DataFrame.

    Parameters:
    -------
    adata       -- An AnnData object which stores the gene expression data and metadata.
    gtf_info    -- A DataFrame containing the mapping from one set of feature names to another.
    from_col    -- A string specifying the column in gtf_info to be mapped from. Default is 'GENEID'.
    to_col      -- A string specifying the column in gtf_info to be mapped to. Default is 'GeneSymbol'.
    inplace     -- A boolean specifying whether to perform the conversion inplace or return a new AnnData object. Default is True.
    verbose     -- A boolean specifying whether to print progress information. Default is True.

    Returns:
    -------
    If inplace is False, returns a new AnnData object with converted feature names.
    """

    if not inplace:
        adata = adata.copy()

    # Filter gtf_info to keep only the gene names found in the anndata object
    if from_col not in gtf_info.columns:
        raise ValueError(f"Column {from_col} not found in gtf_info")
    else:
        gtf_info_filtered = gtf_info[gtf_info[from_col].isin(adata.var_names)]

    if verbose:
        num_found = len(gtf_info_filtered)
        num_total = len(adata.var_names)
        # fraction_found = num_found / num_total
        print(
            f"Fraction of adata.var_names found in gtf_info[{from_col}]: {num_found} out of {num_total}"
        )

    gene_name_mapping = dict(zip(gtf_info[from_col], gtf_info[to_col]))

    adata.var[from_col] = adata.var_names
    adata.var[to_col] = adata.var[from_col].map(gene_name_mapping)

    adata.var.dropna(subset=[to_col], inplace=True)
    adata.var.reset_index(drop=True, inplace=True)

    # mask = adata.var_names.isin(adata.var[to_col].values)
    # adata = adata[:, mask].copy()
    adata.var_names = adata.var[to_col]
    adata.var_names_make_unique()

    if not inplace:
        return adata


# Remove cells with fewer than K neighbors within a distance D
def spatial_singlet_filter(
    adata: ad.AnnData,
    basis="spatial",
    D: float = 100,
    K: int = 10,
    iters: int = 1,
    verbose: bool = True,
    draw_histograms: bool = False,
    inplace: bool = True,
    mode: str = "fast",
):
    """
    This function calculates the Euclidean distances between spatial coordinates,
    finds the number of spatial neighbors within a given distance (D), and filters out cells
    with fewer than a given number of neighbors (K).

    Parameters
    ----------
    adata: sc.AnnData
        An AnnData object. The function expects that this object contains spatial coordinates in adata.obsm['spatial'].

    D: int, default=100
        The maximum Euclidean distance to consider when determining spatial neighbors.

    K: int, default=10
        The minimum number of spatial neighbors required to keep a cell in the output AnnData object.

    iters: int, default=1
        The number of iterations to perform. Each iteration will update the AnnData object with the number of
        spatial neighbors within distance D for each cell. The first iteration will use the original spatial
        coordinates, and subsequent iterations will use the spatial coordinates from the previous iteration.

    verbose: bool, default=True
        If True, print progress messages.

    draw_histograms: bool, default=True
        If True, print histograms showing KNN distributions.

    inplace: bool, default=True
        If True, perform the filtering in-place and modify the input adata object. If False, create a copy of adata and
        perform the filtering on the copy, leaving the original adata object unchanged.

    mode: str, default="fast"
        Mode to choose the implementation. "fast" uses the original implementation that calculates the full distance
        matrix, while "lowmem" uses an iterative approach that calculates neighbors one at a time to save memory.

    Returns
    -------
    adata: ad.AnnData
        An updated AnnData object which only contains cells with more than K neighbors within a distance D.
        Also, this object contains a new field in .obs: "spatial_neighbors_{D}", which contains
        the number of spatial neighbors for each cell within a distance D.
    """
    import scipy.spatial as scisp
    import numpy as np

    if not inplace:
        adata = adata.copy()

    if verbose:
        print(f"{adata.shape[0]} total cells...")

    for i in range(iters):
        ncells = adata.shape[0]
        neighbor_key = f"{basis}_neighbors_{D}"

        if mode == "lowmem":
            # Create a new key in .obs for spatial neighbors within distance D
            neighbor_list = [0] * ncells

            for cell_idx in range(ncells):
                # Get the spatial coordinates of the current cell
                current_coords = adata.obsm[basis][cell_idx]

                # Calculate Euclidean distances to all other cells
                distances = np.linalg.norm(adata.obsm[basis] - current_coords, axis=1)

                # Count the number of neighbors within distance D
                num_neighbors = np.sum(distances < D)

                # Update the number of spatial neighbors in the AnnData object
                # adata.obs.at[adata.obs_names[cell_idx], neighbor_key] = num_neighbors
                neighbor_list[cell_idx] = num_neighbors

            # Filter out cells with fewer than K spatial neighbors within distance D
            adata.obs[neighbor_key] = neighbor_list

            if draw_histograms:
                plt.figure(figsize=(2, 2))
                plt.hist(adata.obs[neighbor_key], bins=100, width=0.2)
                plt.title(f"{i+1} / {iters}")
                plt.xlabel("N_neighbors", fontsize=4)
                plt.ylabel("N_bins", fontsize=4)
                plt.xticks(fontsize=4)
                plt.yticks(fontsize=4)
                plt.axvline(x=K, color="red")
                plt.show()

            adata = adata[adata.obs[neighbor_key] > K, :]

        elif mode == "fast":
            # Calculate Euclidean distances
            basis_distances = scisp.distance.squareform(
                scisp.distance.pdist(adata.obsm[basis])
            )

            # Calculate number of spatial neighbors within distance D for each cell
            adata.obs[neighbor_key] = np.sum(basis_distances < D, axis=0)

            if draw_histograms:
                plt.figure(figsize=(2, 2))
                plt.hist(adata.obs[neighbor_key], bins=100, width=0.2)
                plt.title(f"{i+1} / {iters}")
                plt.xlabel("N_neighbors", fontsize=4)
                plt.ylabel("N_bins", fontsize=4)
                plt.xticks(fontsize=4)
                plt.yticks(fontsize=4)
                plt.axvline(x=K, color="red")
                plt.show()

            # Filter out cells with fewer than K spatial neighbors within distance D
            adata = adata[adata.obs[neighbor_key] > K, :]
        else:
            raise ValueError("Invalid mode. Mode must be 'fast' or 'lowmem'.")

        # Print statements
        if verbose:
            print(f"{i+1} / {iters}: removed {ncells - adata.shape[0]} cells")

        if ncells - adata.shape[0] == 0:
            if verbose:
                print(f"   Finished after {i+1} iterations")
            break

    if not inplace:
        return adata


# Function to segment tissues based on spatial information
def segment_tissues(
    adata, threshold="auto", num_tissues=None, inplace=True, verbose=True
):
    """
    Segment tissues based on spatial information.

    Parameters:
        adata (Anndata): Annotated data object containing spatial coordinates.
        threshold (float or str): Threshold distance for tissue segmentation. If 'auto', the threshold will be calculated based on the expected number of tissues.
        num_tissues (int): Expected number of tissues. Required if threshold='auto'.
        inplace (bool): If True, the tissue labels will be added as adata.obs[f"segment_{threshold}"]. If False, the function will return a modified copy of the input Anndata object.
        verbose (bool): If True, print progress and summary messages.

    Returns:
        Anndata: Annotated data object with tissue labels added as adata.obs[f"segment_{threshold}"].

    """
    if not inplace:
        adata = adata.copy()

    spatial_coords = adata.obsm["spatial"]
    num_cells = len(spatial_coords)

    # Calculate the average distance between cells
    avg_distance = np.mean(
        np.linalg.norm(spatial_coords - np.mean(spatial_coords, axis=0), axis=1)
    )

    # Calculate the threshold based on the expected number of tissues if 'auto' is selected
    if threshold == "auto":
        if num_tissues is None:
            raise ValueError("num_tissues must be provided when threshold='auto'.")
        threshold = num_tissues * avg_distance

    tissue_labels = np.arange(num_cells)
    for i in range(num_cells):
        if verbose and i % (num_cells // 4) == 0:
            print(f"Processing cell {i+1}/{num_cells}")

        distances = np.linalg.norm(spatial_coords - spatial_coords[i], axis=1)
        similar_cells = np.where(distances <= threshold)[0]
        tissue_labels[similar_cells] = tissue_labels[i]

    unique_labels, label_counts = np.unique(tissue_labels, return_counts=True)
    label_mapping = {label: str(i) for i, label in enumerate(unique_labels)}
    tissue_labels = pd.Categorical(tissue_labels.astype(str))
    tissue_labels = tissue_labels.rename_categories(label_mapping)
    adata.obs[f"segment_{threshold}"] = tissue_labels

    if verbose:
        unique_labels, label_counts = np.unique(tissue_labels, return_counts=True)
        print(f"Segmentation completed using threshold={threshold}")
        print(f"Number of identified tissues: {len(unique_labels)}")
        print("Tissue labels summary:")
        for label, count in zip(unique_labels, label_counts):
            print(f"Tissue {label}: {count} cells")

    if not inplace:
        return adata


# Function to add biotype % values to AnnData object
def add_biotypes_pct(
    adata: ad.AnnData,
    biomart: Union[None, pd.DataFrame] = None,  # DataFrame containing gene biotypes
    gene_colname: str = "GeneSymbol",
    biotype_colname: str = "Biotype",
    add_as: str = "obs",  # how percent features should be added
    prefix: str = "pct.",
    scale: int = 100,
    verbose: bool = True,
) -> ad.AnnData:
    """
    This function adds gene biotype percentage values to an AnnData object.

    Args:
        adata (AnnData): The AnnData object containing gene expression data.
        biomart (pd.DataFrame, optional): A DataFrame containing gene biotypes.
        gene_colname (str, optional): Column name in biomart DataFrame for gene identifiers. Default is "GeneSymbol".
        biotype_colname (str, optional): Column name in biomart DataFrame for biotype. Default is "Biotype".
        add_as (str, optional): Determines how percent features should be added. Default is "obs".
        prefix (str, optional): Prefix for column names added to the AnnData object. Default is "pct.".
        scale (int, optional): Determines the scaling for the percentage. Default is 100.
        verbose (bool, optional): Determines whether to print messages during function execution. Default is True.

    Returns:
        AnnData: The original AnnData object with added gene biotype percentage values.
    """

    if biomart is None:
        if verbose:
            print("Need a list of gene biotypes! Nothing done.")
        return adata

    if add_as == "var":
        if verbose:
            print("add_as='var' is not yet implemented")
        return adata

    if verbose:
        print(f"Adding gene biotype percentage values as {add_as} ...")

    biotypes = biomart[biotype_colname].unique()

    for biotype in biotypes:
        # Subset out current gene biotype
        tmp_mart = biomart[biomart[biotype_colname] == biotype]

        # Get unique gene names which are present in the adata object
        tmp_feat = tmp_mart[gene_colname][
            tmp_mart[gene_colname].isin(adata.var_names)
        ].unique()

        if len(tmp_feat) == 0:
            if verbose:
                print(f"  No {biotype} genes found...")
        else:
            if add_as == "obs":
                col_name = prefix + biotype

                # Calculate the percentage
                if issparse(adata.X):
                    gene_pct = adata[:, tmp_feat].X.sum(axis=1) / adata.X.sum(axis=1)
                else:
                    gene_pct = adata[:, tmp_feat].X.sum(axis=1) / adata.X.sum(axis=1)

                # Add the percentage to the AnnData object
                adata.obs[col_name] = np.asarray(gene_pct).flatten()

                # Scale the data if necessary
                if scale == 1:  # [0,1]
                    adata.obs[col_name] /= 100
                elif scale != 100:  # [0,100]
                    if verbose:
                        print("Given scale was not found. Scaling to 100...")
            else:
                if verbose:
                    print("`add_as` option not found... Try again.")

    return adata


# Function to label cells within a region of interest (ROI) polygon
def label_roi_polygon(adata, roi_dict, metadata_col_name="roi", as_string=False):
    """
    Add a column to adata.obs indicating whether each cell is within the defined region of interest (ROI) polygon.

    Parameters:
        adata (AnnData): Anndata object containing spatial coordinates in adata.obsm['spatial'].
        roi_dict (dict): Dictionary containing the vertices of the ROI polygon.
                         Format: {'x': [x1, x2, x3, ...], 'y': [y1, y2, y3, ...]}
        metadata_col_name (str): Name of the new metadata column. Default is 'roi'.
        as_string (bool): If True, store the column values as strings ('True' or 'False').
                          If False, store the column values as booleans (True or False). Default is False.

    Returns:
        None
    """
    roi_path = Path(list(zip(roi_dict["x"], roi_dict["y"])))
    spatial_df = pd.DataFrame(adata.obsm["spatial"], columns=["x", "y"])
    roi_mask = roi_path.contains_points(spatial_df[["x", "y"]].values)

    if as_string:
        roi_mask = pd.Series(roi_mask).map({True: "True", False: "False"})

    adata.obs[metadata_col_name] = roi_mask.values


# Super slow, don't recommend...
def add_tsv_gz_as_layer(
    adata: ad.AnnData,
    tsv_gz_path: str,
    layer_name: str,
    intersect: bool = False,
    transpose: bool = True,
    inplace: bool = True,
):
    """
    Read a tsv.gz file and add it as a new layer to an AnnData object.

    Args:
        adata (AnnData): AnnData object to which the matrix should be added as a new layer.
        tsv_gz_path (str): Path to the tsv.gz file.
        layer_name (str): Name of the layer to be added.
        intersect (bool, optional): Whether to consider only the intersection of observations.
            If False, the union of observations will be used. Default is False.
        transpose (bool, optional): Whether to transpose the matrix before adding it as a layer.
            Default is True.
        inplace (bool, optional): Whether to modify the input AnnData object in place or create a new copy.
            If True, modifications are made to the input AnnData object. If False, a new AnnData object is created
            with the added layer. Default is True.

    Returns:
        None if inplace=True. Returns a new AnnData object with the added layer if inplace=False.
    """
    if not inplace:
        # Create a copy of the AnnData object
        adata = adata.copy()

    # Open the tsv.gz file
    with gzip.open(tsv_gz_path, "rt") as f:
        # Skip the first line
        next(f)

        # Initialize lists to store data
        data_vars = []
        data_obs = []
        data_vals = []

        # Process each line in the file
        for line in f:
            var, obs, val = line.strip().split("\t")
            if obs in adata.obs_names:
                obs_idx = np.where(adata.obs_names == obs)[0][0]
                if float(val) != 0:  # Only add non-zero entries to the sparse matrix
                    data_vars.append(var)
                    data_obs.append(obs_idx)
                    data_vals.append(float(val))

    # Get unique var_names
    var_names = list(set(data_vars))

    # Initialize sparse matrix and transpose the matrix if requested
    if transpose:
        matrix = sp.csr_matrix(
            (data_vals, (data_vars, data_obs)), shape=(adata.n_obs, len(var_names))
        )
    else:
        matrix = sp.csr_matrix(
            (data_vals, (data_obs, data_vars)), shape=(adata.n_obs, len(var_names))
        )

    # Create a new layer in the AnnData object
    adata.layers[layer_name] = matrix

    # Assign column names to the AnnData object
    adata.var_names = var_names

    if not inplace:
        return adata


# Helper fxn to make gene lists unique in the same way anndata does
def make_unique(input_list, delim="-"):
    count_dict = {}
    output_list = []
    for item in input_list:
        count = count_dict.get(item, 0)
        if count > 0:
            unique_item = f"{item}{delim}{count}"
        else:
            unique_item = item
        output_list.append(unique_item)
        count_dict[item] = count + 1
    return output_list


# Load a matrix (saved as an .mtx file) and add as a layer.
def add_mtx_as_layer(
    adata: ad.AnnData,
    mtx_path: str,
    obs_path: str,
    var_path: str,
    layer_name: str,
    feature_column: int = 0,
    transpose: bool = False,
    intersect: bool = False,
    inplace: bool = True,
    verbose: bool = False,
):
    """
    Read an mtx file and add it as a new layer to an AnnData object.

    Args:
        adata (AnnData): AnnData object to which the matrix should be added as a new layer.
        mtx_path (str): Path to the mtx file.
        obs_path (str): Path to row names for the mtx file.
        var_path (str): Path to column names in the mtx file.
        layer_name (str): Name of the layer to be added.
        feature_column (int): Column # for desired feature names (similar to Seurat::ReadMtx)
        transpose (bool, optional): Whether or not to transpose the matrix
        intersect (bool, optional): Whether to consider only the intersection of observations.
            If False, the union of observations will be used. Default is False.
        inplace (bool, optional): Whether to modify the input AnnData object in place or create a new copy.
            If True, modifications are made to the input AnnData object. If False, a new AnnData object is created
            with the added layer. Default is True.

    Returns:
        None if inplace=True. Returns a new AnnData object with the added layer if inplace=False.
    """
    if not inplace:
        # Create a copy of the AnnData object
        adata = adata.copy()

    # Read the matrix file & convert to csr
    # matrix = sp.load_npz(mtx_path)
    with gzip.open(mtx_path, "rb") as file_in:
        matrix = mmread(file_in)
        if transpose:
            matrix = matrix.T
        matrix = matrix.tocsr()

    if verbose:
        print(f"Loaded {matrix.shape[0]} obs and {matrix.shape[1]} vars")

    # Load the associated var_names and row_names from the gzipped .txt files
    # with gzip.open(var_path, 'rt') as var_file:
    #     var_names_mat = var_file.read().splitlines()
    var_names_mat = (
        pd.read_csv(
            var_path, delimiter="\t", dtype=str, header=None, usecols=[feature_column]
        )
        .iloc[:, 0]
        .tolist()
    )

    # Make sure var_names are unique
    var_names_mat = make_unique(var_names_mat)

    with gzip.open(obs_path, "rt") as obs_file:
        obs_names_mat = obs_file.read().splitlines()

    # Filter obs names based on adata obs_names
    # obs_names_found   = list(filter(lambda x: x in adata.obs_names, obs_names_mat))
    # obs_names_missing = list(filter(lambda x: x not in adata.obs_names, obs_names_found))
    obs_names_found = [obs for obs in obs_names_mat if obs in adata.obs_names]
    obs_names_missing = [obs for obs in obs_names_mat if obs not in adata.obs_names]

    # var_names_found   = list(filter(lambda x: x in adata.var_names, var_names_mat))
    # var_names_missing = list(filter(lambda x: x not in adata.var_names, var_names_found))
    var_names_found = [var for var in var_names_mat if var in adata.var_names]
    var_names_missing = [var for var in var_names_mat if var not in adata.var_names]

    if verbose:
        print(
            f"{len(obs_names_found)}/{len(adata.obs_names)} obs from adata found, {len(obs_names_missing)} additional found"
        )

    if verbose:
        print(
            f"{len(var_names_found)}/{len(adata.var_names)} vars from adata found, {len(var_names_missing)} additional found"
        )

    # Filter the matrix and column names based on adata, add zeroes for missing **observations**
    obs_index_dict = {name: index for index, name in enumerate(obs_names_mat)}
    obs_found_indices = [
        obs_index_dict[obs] for obs in adata.obs_names if obs in obs_index_dict
    ]

    matrix_obs_found = matrix[obs_found_indices, :]
    if len(obs_names_missing) > 0:
        matrix_obs_missing = sp.csr_matrix((len(obs_names_missing), matrix.shape[1]))
        matrix_obs_out = sp.vstack((matrix_obs_found, matrix_obs_missing))
    else:
        matrix_obs_out = matrix_obs_found

    # Filter the matrix and column names based on adata, add zeroes for missing **features**
    var_index_dict = {name: index for index, name in enumerate(var_names_mat)}
    var_found_indices = [
        var_index_dict[var] for var in adata.var_names if var in var_index_dict
    ]

    matrix_var_found = matrix_obs_out[:, var_found_indices]
    if len(var_names_missing) > 0:
        matrix_var_missing = sp.csr_matrix(
            (matrix_obs_out.shape[0], len(var_names_missing))
        )
        matrix_out = sp.hstack((matrix_var_found, matrix_var_missing))
    else:
        matrix_out = matrix_var_found

    # Add missing obs_names and var_names to the end of the current lists
    obs_names_found += obs_names_missing
    var_names_found += var_names_missing

    # Now the lengths should match, and you can safely compare the arrays
    if not all(elem in adata.obs_names for elem in obs_names_found) or not all(
        elem in adata.var_names for elem in var_names_found
    ):
        # if any(adata.obs_names != obs_names_found) or any(adata.var_names != var_names_found):
        if verbose:
            print(
                f"Observation names or feature names don't match, getting the order..."
            )
        obs_index_dict = {name: index for index, name in enumerate(adata.obs_names)}
        obs_indices = [obs_index_dict[obs] for obs in adata.obs_names]

        var_index_dict = {name: index for index, name in enumerate(adata.var_names)}
        var_indices = [var_index_dict[var] for var in adata.var_names]

        # Debugging statements
        # print('Shape of matrix_out:', matrix_out.shape)
        # print('Shape of adata:', adata.shape)
        # print('Length of obs_names_found:', len(obs_names_found))
        # print('Length of var_names_found:', len(var_names_found))
        # print('Length of obs_indices:', len(obs_indices))
        # print('Length of var_indices:', len(var_indices))

        matrix_out = matrix_out[np.ix_(obs_indices, var_indices)]

    # Create a new layer in the AnnData object
    adata.layers[layer_name] = matrix_out

    if not inplace:
        return adata


# Rotate the spatial coordinates in an AnnData object
def rotate_embedding(adata, basis, theta, inplace=True):
    """
    Rotate spatial coordinates stored in an anndata object for a specific embedding basis clockwise by a given angle.

    Parameters:
    adata (anndata.AnnData): Anndata object containing spatial coordinates in adata.obsm['spatial'].
    basis (str): Name of the embedding basis in which the spatial coordinates are stored.
    theta (float): Angle of rotation in degrees.
    inplace (bool): If True, the rotation is performed in place and modifies the original anndata object.
                    If False (default), a new anndata object with rotated spatial coordinates is returned.

    Returns:
    anndata.AnnData or None: If inplace=True, returns None. If inplace=False, returns a new anndata object with
                             rotated spatial coordinates.
    """

    if not inplace:
        adata = adata.copy()

    if basis not in adata.obsm_keys():
        raise ValueError(
            f"The specified embedding basis '{basis}' is not found in adata.obsm."
        )

    # Convert theta from degrees to radians
    theta_rad = np.radians(theta)

    # Extract spatial coordinates from adata
    spatial_coords = adata.obsm[basis]

    # Get the center of the spatial coordinates
    center = np.mean(spatial_coords, axis=0)

    # Define the 2D rotation matrix
    rotation_matrix = np.array(
        [
            [np.cos(theta_rad), np.sin(theta_rad)],
            [-np.sin(theta_rad), np.cos(theta_rad)],
        ]
    )

    # Apply rotation to each spatial coordinate relative to the center
    rotated_coords = np.dot(spatial_coords - center, rotation_matrix) + center

    # Update the original anndata object with rotated spatial coordinates
    adata.obsm[basis] = rotated_coords

    if not inplace:
        return adata


# Rescale spatial coordinates in an AnnData object
def rescale_embedding(
    adata, basis, scale_factor, move_to_origin=False, inplace=True, verbose=False
):
    """
    Rescale spatial coordinates in an anndata object.

    Parameters:
        adata (anndata.AnnData): Annotated data object.
        basis (str): Key for the basis containing the spatial coordinates to rescale.
        scale_factor (float): Scale factor to apply to the spatial coordinates.
        move_to_origin (bool, optional): If true, translocate coordinates so that the min() of each axis is zero
        inplace (bool, optional): If True, perform the operation in place.
                                  If False, return a new annotated data object with the rescaled spatial coordinates.
                                  Default is True.
        verbose (bool, optional): If True, print updates.

    Returns:
        anndata.AnnData or None: If inplace is True, returns None.
                                 If inplace is False, returns a new annotated data object with the rescaled spatial coordinates.
    """
    if not inplace:
        adata = adata.copy()

    if basis not in adata.obsm_keys():
        raise ValueError(
            f"The basis '{basis}' is not present in the 'obsm' of the anndata object."
        )

    # Get the spatial coordinates from the 'obsm' attribute
    coords = adata.obsm[basis]

    # Handle iso/anisotropic scaling
    if type(scale_factor) in [float, int]:
        scale_factor = [scale_factor] * coords.shape[1]

    # Translocate the coordinates so that the leftmost
    if move_to_origin:
        if verbose:
            print(f"Moving coordinates to the origin")
        for i in range(len(scale_factor)):
            coords[:, i] = coords[:, i] - min(coords[:, i])

    # Rescale the spatial coordinatesif verbose:
    if verbose:
        print(
            f"   x range: {min(adata.obsm[basis][:,0])}, {max(adata.obsm[basis][:,0])}"
        )
        print(
            f"   y range: {min(adata.obsm[basis][:,1])}, {max(adata.obsm[basis][:,1])}"
        )

    if len(scale_factor) == 1:  # isotropic scaling
        coords = coords * scale_factor
    elif len(scale_factor) > 1:  # anisotropic scaling
        if verbose:
            print(f"Performing anisotropic scaling...")

        for i in range(len(scale_factor)):
            if verbose:
                print(f"     Using scale factor of {scale_factor[i]} for dimension {i}")
            coords[:, i] = coords[:, i] * scale_factor[i]

    # Update the 'obsm' attribute with the rescaled spatial coordinates
    adata.obsm[basis] = coords

    # Final verbose print
    if verbose:
        print(f"Finished scaling:")
        print(
            f"   x range: {min(adata.obsm[basis][:,0])}, {max(adata.obsm[basis][:,0])}"
        )
        print(
            f"   y range: {min(adata.obsm[basis][:,1])}, {max(adata.obsm[basis][:,1])}"
        )

    if not inplace:
        return adata


# Get the most abundant genes in an AnnData object
def top_n_genes(adata, n):
    """
    Identify the top N genes with the highest sum of expression values in an AnnData object.

    Parameters:
    adata (AnnData): An AnnData object storing the gene expression data and metadata.
    n (int): The number of top genes to return.

    Returns:
    list: A list of top N genes sorted by their sum of expression values.
    """
    sorted_genes = adata.var_names[np.argsort(adata.X.sum(axis=0))[::-1]]
    return sorted_genes[:n].tolist()

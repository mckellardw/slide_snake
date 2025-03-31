# Utility functions for use with scanpy
import os
import gzip
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np

from matplotlib.path import Path
import matplotlib.pyplot as plt

# from typing import Union

import scipy.sparse as sp
from scipy.sparse import issparse
from scipy.io import mmread


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


def load_gtf_to_dataframe(
    gtf_file, feature_type="all", seqname_filter="all", unwrap_attributes=True
):
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
        for attribute in attributes_str.strip().split(";"):
            if attribute.strip():
                key, value = attribute.strip().split(" ", 1)
                attributes[key] = value.strip('"')
        return attributes

    # Open the file (support for gzipped files)
    open_func = gzip.open if gtf_file.endswith(".gz") else open
    with open_func(gtf_file, "rt") as f:
        # Read the GTF file into a DataFrame
        gtf_data = []
        attributes_data = []
        for line in f:
            if line.startswith("#"):
                continue  # Skip comment lines
            fields = line.strip().split("\t")
            if feature_type != "all" and fields[2] != feature_type:
                continue  # Skip features that don't match the specified type
            if seqname_filter != "all" and fields[0] != seqname_filter:
                continue  # Skip sequence names that don't match the filter
            attributes = parse_attributes(fields[8])
            gtf_data.append(fields[:8])
            attributes_data.append(attributes)

        # Create the main DataFrame
        columns = [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
        ]
        gtf_df = pd.DataFrame(gtf_data, columns=columns)

        if unwrap_attributes:
            # Create a DataFrame for attributes and concatenate with the main DataFrame
            attributes_df = pd.DataFrame(attributes_data)
            gtf_df = pd.concat([gtf_df, attributes_df], axis=1)

    return gtf_df


# Pattern matching for gene names.
## Returns the genes sorted by expression (high to low) by default
def grep_genes(
    adata,
    pattern="",
    layer=None,
    filter_pattern=None,
    sort_by="expression",
    verbose=True,
):
    """
    Search for genes in an AnnData object that match a given pattern.

    Parameters:
    adata (anndata.AnnData): The AnnData object containing the gene expression data.
    pattern (str or list of str, optional): The pattern or patterns to search for in the gene names. Default is "".
    assay (str, optional): The layer to use. Default is None, which uses "counts".
    filter_pattern (str or list of str, optional): A pattern or patterns to filter out from the gene names. Default is None.
    sort_by (str, optional): How to sort the output genes. Can be "expression" (sort by total expression across all cells) or "abc" (sort alphabetically). Default is "expression".
    verbose (bool, optional): Whether to print progress messages. Default is True.

    Returns:
    list of str: The names of the genes that match the pattern and do not match the filter_pattern, sorted according to sort_by.
    """

    if layer is not None:
        print("WARNING: haven't implemented this yet, using `X`")

    genes = adata.var_names

    if verbose:
        print(f"Found {len(genes)} total features...")

    if pattern == "":
        if verbose:
            out_genes = genes
            print("Returning all of the genes!")
    elif isinstance(pattern, list):
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


def grep_var(
    adata: ad.AnnData,
    pattern="",
    var_colname: str = None,
    layer=None,
    filter_pattern=None,
    sort_by="expression",
    verbose=True,
):
    """
    Search for entries in adata.var that match a given pattern.

    Parameters:
    adata (anndata.AnnData): The AnnData object containing the gene expression data.
    pattern (str or list of str, optional): The pattern or patterns to search for. Default is "".
    var_colname (str, optional): The column in adata.var to search. If None, defaults to var_names. Default is None.
    layer (str, optional): The layer to use for expression-based sorting. Default is None, which uses .X.
    filter_pattern (str or list of str, optional): A pattern or patterns to filter out from the results. Default is None.
    sort_by (str, optional): How to sort the output. Can be "expression" (sort by total expression across all cells) or "abc" (sort alphabetically). Default is "expression".
    verbose (bool, optional): Whether to print progress messages. Default is True.

    Returns:
    list of str: The names or values in the specified column that match the pattern and do not match the filter_pattern, sorted according to sort_by.
    """
    if var_colname is None:
        if verbose:
            print("No var_colname provided. Defaulting to var_names.")
        search_data = adata.var_names
    else:
        if var_colname not in adata.var.columns:
            raise ValueError(f"Column '{var_colname}' not found in adata.var.")
        if not pd.api.types.is_categorical_dtype(adata.var[var_colname]) and not pd.api.types.is_string_dtype(adata.var[var_colname]):
            raise ValueError(f"Column '{var_colname}' must be categorical or strings.")
        search_data = adata.var[var_colname]

    if verbose:
        print(f"Searching in column '{var_colname or 'var_names'}'...")

    if pattern == "":
        if verbose:
            print("No pattern provided. Returning all entries.")
        out_entries = search_data
    elif isinstance(pattern, list):
        if verbose:
            print("Looking for multiple patterns...")
        out_entries = [entry for entry in search_data if any(pat in str(entry) for pat in pattern)]
    else:
        if verbose:
            print(f"Looking for pattern '{pattern}'...")
        out_entries = [entry for entry in search_data if pattern in str(entry)]

    if len(out_entries) == 0:
        if verbose:
            print("No matches found.")
        return []

    if filter_pattern is not None:
        if verbose:
            print(f"Filtering out entries containing '{filter_pattern}'...")
        out_entries = [
            entry for entry in out_entries if not any(pat in str(entry) for pat in filter_pattern)
        ]

    if sort_by == "expression" and var_colname is None:
        if verbose:
            print("Sorting by expression...")
        if layer:
            gene_sums = np.asarray(adata[:, out_entries].layers[layer].sum(axis=0)).squeeze()
        else:
            gene_sums = np.asarray(adata[:, out_entries].X.sum(axis=0)).squeeze()
        out_entries = [x for _, x in sorted(zip(gene_sums, out_entries), reverse=True)]
    elif sort_by == "abc":
        if verbose:
            print("Sorting alphabetically...")
        out_entries = sorted(out_entries)

    return out_entries


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
    from_col: str = "GeneID",
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
    from_col    -- A string specifying the column in gtf_info to be mapped from. Default is 'GeneID'.
    to_col      -- A string specifying the column in gtf_info to be mapped to. Default is 'GeneSymbol'.
    inplace     -- A boolean specifying whether to perform the conversion inplace or return a new AnnData object. Default is True.
    verbose     -- A boolean specifying whether to print progress information. Default is True.

    Returns:
    -------
    If inplace is False, returns a new AnnData object with converted feature names.
    """

    if not inplace:
        adata = adata.copy()

    original_shape = adata.shape

    if verbose:
        print(f"Original AnnData shape: {original_shape}")
        print(f"Original var_names length: {len(adata.var_names)}")
        print(f"gtf_info shape: {gtf_info.shape}")

    # Check if required columns are in gtf_info
    if from_col not in gtf_info.columns or to_col not in gtf_info.columns:
        raise ValueError(f"Columns {from_col} and/or {to_col} not found in gtf_info")

    # Create a mapping dictionary
    gene_name_mapping = dict(zip(gtf_info[from_col], gtf_info[to_col]))

    # Create a new column with mapped names, keeping original names for unmapped genes
    adata.var["original_names"] = adata.var_names
    adata.var["new_names"] = adata.var_names.map(lambda x: gene_name_mapping.get(x, x))

    # Count and report on mapping results
    if verbose:
        total_genes = len(adata.var_names)
        mapped_genes = sum(adata.var["new_names"] != adata.var["original_names"])
        print(f"Total genes: {total_genes}")
        print(f"Mapped genes: {mapped_genes}")
        print(f"Unmapped genes: {total_genes - mapped_genes}")

    # Update var_names
    adata.var_names = adata.var["new_names"]

    # Make var_names unique
    adata.var_names_make_unique()

    # Clean up temporary columns
    adata.var.drop(columns=["original_names", "new_names"], inplace=True)

    if verbose:
        print(f"Final AnnData shape: {adata.shape}")
        print(f"Final var_names length: {len(adata.var_names)}")
        if adata.shape != original_shape:
            print("WARNING: AnnData shape has changed!")
            print(
                "Genes lost:", set(adata.var["original_names"]) - set(adata.var_names)
            )

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
    biomart: pd.DataFrame,  # DataFrame containing gene biotypes
    gene_colname: str = "GeneSymbol",
    biotype_colname: str = "Biotype",
    add_as: str = "obs",  # how percent features should be added
    prefix: str = "pct_",
    scale: int = 100,
    verbose: bool = True,
    layer: str = None,  # Specify the layer to use
    inplace: bool = True,  # Whether to modify the input AnnData object in place
) -> ad.AnnData:
    """
    This function adds gene biotype percentage values to an AnnData object.

    Args:
        adata (AnnData): The AnnData object containing gene expression data.
        biomart (pd.DataFrame, optional): A DataFrame containing gene biotypes.
        gene_colname (str, optional): Column name in biomart DataFrame for gene identifiers. Default is "GeneSymbol".
        biotype_colname (str, optional): Column name in biomart DataFrame for biotype. Default is "Biotype".
        add_as (str, optional): Determines how percent features should be added. Default is "obs".
        prefix (str, optional): Prefix for column names added to the AnnData object. Default is "pct_".
        scale (int, optional): Determines the scaling for the percentage. Default is 100.
        verbose (bool, optional): Determines whether to print messages during function execution. Default is True.
        layer (str, optional): The layer in the AnnData object to use for computation. Default is None (uses .X).
        inplace (bool, optional): Whether to modify the input AnnData object in place. Default is True.

    Returns:
        AnnData: The modified AnnData object if inplace=False, otherwise None.
    """
    if not inplace:
        adata = adata.copy()

    if biomart is None:
        if verbose:
            print("Need a list of gene biotypes! Nothing done.")
        return adata if not inplace else None

    if add_as == "var":
        if verbose:
            print("add_as='var' is not yet implemented")
        return adata if not inplace else None

    # Check if the specified layer exists
    if layer is not None and layer not in adata.layers:
        if verbose:
            print(
                f"The specified layer '{layer}' does not exist in the AnnData object. Using .X instead."
            )
        layer = None

    if verbose:
        print(f"Adding gene biotype percentage values as {add_as} ...")
        if layer:
            print(f"Using layer: {layer}")
        else:
            print("Using default expression matrix (.X)")

    biotypes = list(biomart[biotype_colname].unique())
    if verbose:
        print(f"Found {len(biotypes)} biotypes..")

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

                # Calculate the percentage based on the specified layer or .X
                if layer:
                    gene_pct = adata[:, tmp_feat].layers[layer].sum(
                        axis=1
                    ) / adata.layers[layer].sum(axis=1)
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

    return adata if not inplace else None


def add_var_pct(
    adata: ad.AnnData,
    var_colname: str = "gene_type",
    add_as: str = "obs",  # how percent features should be added
    prefix: str = "pct_",
    scale: int = 100,
    verbose: bool = True,
    layer: str = None,  # Specify the layer to use
    inplace: bool = True,  # Whether to modify the input AnnData object in place
) -> ad.AnnData:
    """
    This function computes percentages based on a column in adata.var and adds them to the AnnData object.

    Args:
        adata (AnnData): The AnnData object containing gene expression data.
        var_colname (str, optional): Column name in adata.var to group features by. Default is "gene_type".
        add_as (str, optional): Determines how percent features should be added. Default is "obs".
        prefix (str, optional): Prefix for column names added to the AnnData object. Default is "pct_".
        scale (int, optional): Determines the scaling for the percentage. Default is 100.
        verbose (bool, optional): Determines whether to print messages during function execution. Default is True.
        layer (str, optional): The layer in the AnnData object to use for computation. Default is None (uses .X).
        inplace (bool, optional): Whether to modify the input AnnData object in place. Default is True.

    Returns:
        AnnData: The modified AnnData object if inplace=False, otherwise None.
    """
    if not inplace:
        adata = adata.copy()

    if var_colname not in adata.var.columns:
        if verbose:
            print(f"Column '{var_colname}' not found in adata.var. Nothing done.")
        return adata if not inplace else None

    # Check if the specified layer exists
    if layer is not None and layer not in adata.layers:
        if verbose:
            print(
                f"The specified layer '{layer}' does not exist in the AnnData object. Using .X instead."
            )
        layer = None

    if verbose:
        print(f"Adding percentages based on '{var_colname}' as {add_as} ...")
        if layer:
            print(f"Using layer: {layer}")
        else:
            print("Using default expression matrix (.X)")

    categories = adata.var[var_colname].unique()
    if verbose:
        print(f"Found {len(categories)} unique categories in '{var_colname}'.")

    for category in categories:
        # Subset out current category
        tmp_feat = adata.var_names[adata.var[var_colname] == category]

        if len(tmp_feat) == 0:
            if verbose:
                print(f"  No features found for category '{category}'...")
        else:
            if add_as == "obs":
                col_name = prefix + category

                # Calculate the percentage based on the specified layer or .X
                if layer:
                    total_counts = adata.layers[layer].sum(axis=1)
                    gene_pct = adata[:, tmp_feat].layers[layer].sum(axis=1) / np.where(
                        total_counts == 0, np.nan, total_counts
                    )
                else:
                    total_counts = adata.X.sum(axis=1)
                    gene_pct = adata[:, tmp_feat].X.sum(axis=1) / np.where(
                        total_counts == 0, np.nan, total_counts
                    )

                # Replace NaN values with 0
                gene_pct = np.nan_to_num(gene_pct)

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

    return adata if not inplace else None


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
    assume_order: bool = False,
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
        feature_column (int): Column # for desired feature names (similar to Seurat::ReadMtx).
        transpose (bool, optional): Whether or not to transpose the matrix.
        intersect (bool, optional): Whether to consider only the intersection of observations and features.
            If False, the union of observations and features will be used. Default is False.
        assume_order (bool, optional): If True, assumes the adata object and the new matrix have features in the same order.
            If True, directly loads the matrix without reordering. Default is True.
        inplace (bool, optional): Whether to modify the input AnnData object in place or create a new copy.
            If True, modifications are made to the input AnnData object. If False, a new AnnData object is created
            with the added layer. Default is True.
        verbose (bool, optional): Whether to print additional information. Default is False.

    Returns:
        None if inplace=True. Returns a new AnnData object with the added layer if inplace=False.
    """
    if not inplace:
        adata = adata.copy()

    with gzip.open(mtx_path, "rb") as file_in:
        matrix = mmread(file_in)
        if transpose:
            matrix = matrix.T
        matrix = matrix.tocsr()

    if verbose:
        print(f"Loaded {matrix.shape[0]:,} obs and {matrix.shape[1]:,} vars")

    if assume_order:
        if verbose:
            print("Assuming the same feature order between adata and the new matrix.")
        adata.layers[layer_name] = matrix
        return adata if not inplace else None

    var_names_mat = (
        pd.read_csv(
            var_path, delimiter="\t", dtype=str, header=None, usecols=[feature_column]
        )
        .iloc[:, 0]
        .tolist()
    )
    var_names_mat = make_unique(var_names_mat)

    with gzip.open(obs_path, "rt") as obs_file:
        obs_names_mat = obs_file.read().splitlines()

    if intersect:
        obs_names_common = list(set(adata.obs_names) & set(obs_names_mat))
        var_names_common = list(set(adata.var_names) & set(var_names_mat))

        if verbose:
            print(f"Subsetting to {len(obs_names_common):,} obs and {len(var_names_common):,} vars")

        # Subset and re-order the AnnData object and the matrix
        adata = adata[obs_names_common, var_names_common]
        obs_indices_mat = [obs_names_mat.index(obs) for obs in obs_names_common]
        var_indices_mat = [var_names_mat.index(var) for var in var_names_common]
        matrix = matrix[obs_indices_mat, :][:, var_indices_mat]
    else:
        obs_names_out = list(set(adata.obs_names) | set(obs_names_mat))
        var_names_out = list(set(adata.var_names) | set(var_names_mat))

        if verbose:
            print(f"Final matrix shape: {len(obs_names_out):,} obs and {len(var_names_out):,} vars")

        # Create a new matrix with the correct dimensions
        new_matrix = sp.csr_matrix(
            (len(obs_names_out), len(var_names_out)), dtype=matrix.dtype
        )

        # Fill in the values from the original matrix
        obs_index_dict = {name: index for index, name in enumerate(obs_names_mat)}
        var_index_dict = {name: index for index, name in enumerate(var_names_mat)}

        for i, obs in enumerate(obs_names_out):
            if obs in obs_index_dict:
                orig_obs_idx = obs_index_dict[obs]
                for j, var in enumerate(var_names_out):
                    if var in var_index_dict:
                        orig_var_idx = var_index_dict[var]
                        new_matrix[i, j] = matrix[orig_obs_idx, orig_var_idx]

        # Update adata object
        missing_obs = list(set(obs_names_out) - set(adata.obs_names))
        missing_vars = list(set(var_names_out) - set(adata.var_names))

        if missing_obs:
            adata = adata.concatenate(
                ad.AnnData(
                    X=sp.csr_matrix((len(missing_obs), adata.n_vars)),
                    obs=pd.DataFrame(index=missing_obs),
                ),
                join="outer",
            )

        if missing_vars:
            adata = adata.concatenate(
                ad.AnnData(
                    X=sp.csr_matrix((adata.n_obs, len(missing_vars))),
                    var=pd.DataFrame(index=missing_vars),
                ),
                join="outer",
                axis=1,
            )

        # Ensure the order of observations and variables matches
        adata = adata[obs_names_out, var_names_out]

        # Assign the new matrix
        matrix = new_matrix

    # Add the new layer
    adata.layers[layer_name] = matrix

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


def add_barcode_tally_to_anndata(
    tsv_path, adata, metadata_key="barcode_count", add_missing=False
):
    """
    Reads a TSV file containing read_IDs and barcodes, tallies the barcodes,
    and adds the counts as metadata to an AnnData object.

    Parameters:
    -----------
    tsv_path : str
        Path to the TSV file containing read_IDs and barcodes
    adata : AnnData
        AnnData object to add the metadata to
    metadata_key : str, optional
        Key to use for the metadata in adata.obs (default: 'barcode_count')
    add_missing : bool, optional
        If True, adds missing barcodes to the AnnData object
        If False, skips barcodes that aren't in the AnnData object (default: False)

    Returns:
    --------
    None (modifies adata in place)
    """
    # Read TSV file without header
    df = pd.read_csv(tsv_path, sep="\t", header=None, names=["read_id", "barcode"])

    # Count occurrences of each barcode
    barcode_counts = df["barcode"].value_counts()

    if add_missing:
        # Get unique barcodes from both sources
        tsv_barcodes = set(barcode_counts.index)
        adata_barcodes = set(adata.obs_names)
        all_barcodes = sorted(tsv_barcodes.union(adata_barcodes))

        # Create a new AnnData with all barcodes
        # First, create a DataFrame with zeros for all features
        new_data = pd.DataFrame(0, index=all_barcodes, columns=adata.var_names)

        # Fill in the existing data
        new_data.loc[adata.obs_names] = (
            adata.X.toarray() if issparse(adata.X) else adata.X
        )

        # Create new AnnData object
        new_adata = ad.AnnData(new_data)

        # Copy over existing obs and var annotations
        for col in adata.obs.columns:
            new_adata.obs[col] = pd.Series(index=all_barcodes)
            new_adata.obs.loc[adata.obs_names, col] = adata.obs[col]

        new_adata.var = adata.var.copy()

        # Update the reference to adata
        adata = new_adata

    # Initialize counts for all cells in adata with 0
    adata.obs[metadata_key] = 0

    # Update counts only for barcodes that exist in adata
    existing_barcodes = barcode_counts.index.intersection(adata.obs_names)
    adata.obs.loc[existing_barcodes, metadata_key] = barcode_counts[existing_barcodes]

    # Convert counts to integer type
    adata.obs[metadata_key] = adata.obs[metadata_key].astype(int)

    return adata


#

import PIL

PIL.Image.MAX_IMAGE_PIXELS = 100_000_000

from skimage.transform import rescale, resize
from skimage.io import imread, imshow, imsave


def rescale_annData(adata, adata2):
    """
    Rescale the spatial coordinates and image of an AnnData object to match another AnnData object.

    Parameters:
    adata (AnnData): The AnnData object to be rescaled.
    adata2 (AnnData): The AnnData object to match the scale.

    Returns:
    AnnData: A new AnnData object with rescaled spatial coordinates and image.
    """
    adata1 = adata.copy()
    img1 = adata1.uns["spatial"]["HE"]
    p_ma_y, p_ma_x = img1.shape[1], img1.shape[0]
    img2 = adata2.uns["spatial"]["HE"]
    n_ma_y, n_ma_x = img2.shape[1], img2.shape[0]
    if len(adata1.obsm["spatial"].T) > 2:
        new_spatial = adata1.obsm["spatial"].copy()
    else:
        new_spatial = adata1.obsm["spatial"].T.copy()
    new_spatial[0] = n_ma_y / p_ma_y * new_spatial[0]
    new_spatial[1] = n_ma_x / p_ma_x * new_spatial[1]
    if len(adata1.obsm["spatial"].T) > 2:
        adata1.obsm["spatial"] = new_spatial
    else:
        adata1.obsm["spatial"] = new_spatial.T

    adata1.uns["spatial"]["HE"] = resize(adata1.uns["spatial"]["HE"], (n_ma_x, n_ma_y))
    return adata1


def load_annData(folder, show_img=True):
    """
    Load AnnData object from a folder containing count matrix, barcodes, image, and annotation files.

    Parameters:
    folder (str): Path to the folder containing the files.
    show_img (bool): Whether to display the image. Default is True.

    Returns:
    AnnData: The loaded AnnData object.
    """
    # RAW Data 50um
    count_matrix, barcodes, img, annotation = None, None, None, None
    for file in os.listdir(folder):
        if file.__contains__("count"):
            count_matrix_path = folder + file
            if file.endswith(".tsv"):
                sep = "\t"
            else:
                sep = ","
            count_matrix = pd.read_csv(
                count_matrix_path, index_col=0, header=0, sep=sep
            )
        if file.__contains__("bc"):
            barcodes_path = folder + file
            barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")
            try:
                barcodes.columns = ["barcodes", "x", "y", "z"]
            except ValueError:
                barcodes.columns = ["barcodes", "x", "y"]

        if (
            file.__contains__("png")
            or file.__contains__("tif")
            or file.__contains__("jpg")
        ):
            image_path = folder + file
            img = rescale(imread(image_path), 1 / 2, channel_axis=-1)

        if file.__contains__("annotation"):
            path = folder + file
            annotation = pd.read_csv(path, index_col=0)

    list_col = [
        col for col in count_matrix.columns if col in barcodes["barcodes"].to_numpy()
    ]
    barcodes = barcodes[
        [barcode in count_matrix.columns for barcode in barcodes["barcodes"].to_numpy()]
    ]
    count_matrix = count_matrix[list_col]
    if show_img:
        imshow(img)

    adata = ad.AnnData(count_matrix[barcodes["barcodes"]].T)
    adata.vars = pd.Series(count_matrix.T.columns)

    ma_y, ma_x = img.shape[1], img.shape[0]
    list_ = [
        ma_y * (barcodes["x"]) / np.max(barcodes["x"]),
        ma_x * (50 - barcodes["y"]) / np.max(barcodes["y"]),
    ]
    adata.obsm["spatial"] = np.array(list_).T
    adata.uns["spatial"] = {}
    adata.uns["spatial"]["images"] = {}
    adata.uns["spatial"]["images"]["HE"] = img
    adata.uns["Annotation"] = annotation
    return adata


def load_visium(folder):
    """
    Load Visium data from a folder containing raw feature matrix, image, and annotation files.

    Parameters:
    folder (str): Path to the folder containing the files.

    Returns:
    AnnData: The loaded Visium AnnData object.
    """
    visium, bc, visium_img, annotation = None, None, None, None
    for file in os.listdir(folder):
        if file.__contains__("raw_feature_bc_matrix") and file.endswith(".h5"):
            h5_path = folder + file
            visium = sc.read_10x_h5(h5_path)
        if file.__contains__("png") or file.__contains__("tif"):
            image_path = folder + file
            visium_img = rescale(imread(image_path), 1 / 4, channel_axis=-1)
        if file.__contains__("annotation"):
            path = folder + file
            annotation = pd.read_csv(path, index_col=0)
        if file.__contains__("spatial"):
            try:
                for file2 in os.listdir(folder + file):
                    if file2.endswith("json"):
                        json_file_path = folder + file + "/" + file2
                    if file2.__contains__("tissue_positions"):
                        barcodes_path = folder + file + "/" + file2
                        bc = pd.read_csv(barcodes_path, index_col=0)
            except:
                None

    imshow(visium_img)
    visium.uns["spatial"] = {}
    visium.uns["spatial"]["HE"] = visium_img
    visium.X = np.array(visium.X.todense())
    visium = visium[bc.index].copy()
    visium.obsm["spatial"] = bc
    visium.obsm["spatial"].columns = [*visium.obsm["spatial"].columns[:-2], 1, 0]
    visium.obsm["spatial"][0] = visium.obsm["spatial"][0] / 4
    visium.obsm["spatial"][1] = visium.obsm["spatial"][1] / 4
    visium.uns["Annotation"] = annotation
    return visium


# #To Celery (And beyond)
def to_celery(adata, out_folder):
    """
    Export AnnData object to a format compatible with Celery.

    Parameters:
    adata (AnnData): The AnnData object to be exported.
    out_folder (str): Path to the output folder.

    Returns:
    None
    """

    def compose_alpha(image_with_alpha):
        """
        Compose an image with an alpha channel into an image without an alpha channel.

        Parameters:
        image_with_alpha (ndarray): The image with an alpha channel.

        Returns:
        ndarray: The composed image without an alpha channel.
        """
        image_with_alpha = image_with_alpha.astype(np.float32)
        image, alpha = image_with_alpha[..., :3], image_with_alpha[..., 3:] / 255.0
        image = image * alpha + (1.0 - alpha) * 255.0
        image = image.astype(np.uint8)
        return image

    if len(adata.uns["spatial"]["HE"]).shape == 4:
        img = compose_alpha(adata.uns["spatial"]["HE"])
    else:
        img = adata.uns["spatial"]["HE"]
    np.save(
        out_folder + "/coord_inside_grid.npy",
        np.array([adata.obsm["spatial"][1], adata.obsm["spatial"][0]]).T,
    )
    imsave(out_folder + "/aligned.png", img)
    pd.DataFrame(data=adata.X, index=adata.var, columns=adata.obs).to_csv(
        out_folder + "/sequencing_inside.csv"
    )
    pd.DataFrame(adata.X.var.index.to_numpy()).to_csv(
        out_folder + "/gene_names.csv", index=False, header=False
    )

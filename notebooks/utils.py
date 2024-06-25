import anndata as ad
import os
import pandas as pd
from skimage.transform import rescale, resize
from skimage.io import imread, imshow
import numpy as np
import PIL
from skimage.io import imsave
import json
import h5py
import numpy as np
import scanpy as sc
import matplotlib as mpl

PIL.Image.MAX_IMAGE_PIXELS = 100_000_000


def rescale_annData(adata, adata2):
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


def processing(adata, inplace=True):
    min_counts = 100
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.normalize_total(adata, target_sum=1e4, inplace=inplace)
    print("Normalized")
    sc.pp.log1p(adata)
    print("log1pd")
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=1000)
    print("highly variable gened")
    sc.pp.pca(adata)
    print("pcad")
    sc.pp.neighbors(adata)
    print("neighbored")
    sc.tl.leiden(adata, key_added="leiden")
    print("leidend")
    return adata


# #To Celery (And beyond)
def to_celery(adata, out_folder):
    def compose_alpha(image_with_alpha):

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


# Plotting functions for use with scanpy
import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from anndata import AnnData

# import anndata as ad


# Knee plot to quality check UMI counts for single-cell data
def knee_plot(
    ADATA,
    x_lim=[0, 20000],
    line_width=2,
    line_color="b",
    title="Knee plot",
    verbose=False,
):
    import matplotlib.pyplot as plt

    expected_num_cells = 10000

    knee = np.sort((np.array(ADATA.X.sum(axis=1))).flatten())[::-1]

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.loglog(knee, range(len(knee)), linewidth=line_width, color=line_color)
    #     ax.axvline(x=knee[expected_num_cells], linewidth=3, color="k")
    #     ax.axhline(y=expected_num_cells, linewidth=3, color="k")

    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Set of Barcodes")
    ax.set_title(title)

    plt.xlim(x_lim)
    plt.grid(True, which="both")
    plt.show()


# Faceted plot for any embedding
# scanpy github issue reference- https://github.com/scverse/scanpy/issues/955
def facet_embedding(
    adata, clust_key, basis, size=60, frameon=False, legend_loc=None, **kwargs
):
    # import scanpy as sc

    tmp = adata.copy()

    for i, clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype("category")
        tmp.uns[clust + "_colors"] = ["#d3d3d3", adata.uns[clust_key + "_colors"][i]]

    sc.pl.embedding(
        tmp,
        groups=tmp.obs[clust].cat.categories[1:].values,
        color=adata.obs[clust_key].cat.categories.tolist(),
        basis=basis,
        size=size,
        frameon=frameon,
        legend_loc=legend_loc,
        **kwargs,
    )


##David's code
# scanpy version of seuListPlot - grids of plots for a single embedding
from mpl_toolkits.axes_grid1 import make_axes_locatable


# import mpl_toolkits.axes_grid1.inset_locator as mpl_il
def plot_grid_of_embeddings(
    adata_dict_,
    color,
    ncols=None,
    titles=None,
    figsize=(10, 10),
    na_color="lightgrey",
    colorbar_width=10,
    colorbar_pad=0.05,
    same_scale=False,
    rescale_images=False,
    filename=None,
    **kwargs,
):
    """
    Plot a grid of embeddings for multiple Anndata objects.
    Parameters:
        adata_dict (dict): Dictionary of Anndata objects. The keys represent the plot titles and the values are the corresponding Anndata objects.
        color (list): List of feature names to use for coloring the embeddings.
        ncols (int): Number of columns in the grid. If None, it is set to the number of entries in adata_dict.
        figsize (tuple): Figure size (width, height).
        same_scale (bool): If True, use the same color scale for plots showing the same feature.
        **kwargs: Additional parameters to be passed to scanpy.pl.embedding.

    """
    adata_dict = adata_dict_.copy()
    # param checks/sets
    if ncols is None:
        ncols = len(adata_dict)

    if titles is None:
        titles = list(adata_dict.keys())

    nrows = len(list(color))

    if len(color) >= 1:
        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=figsize,
            sharex=True,
            sharey="row",
            # gridspec_kw={'hspace': 0}
        )
        axes = axes.flatten()
    else:
        print("Need to specify `color`")
    if rescale_images:
        list_shape = [
            img.shape
            for img in [
                adata.uns["spatial"]["HE"] for (title, adata) in adata_dict.items()
            ]
        ]
        big_img = np.argmax(list_shape)
        for j, (title, adata) in enumerate(adata_dict.items()):
            if j != big_img:
                key = [
                    key
                    for key in adata_dict.keys()
                    if (adata_dict[key].X.shape == adata.X.shape)
                ][0]
                adata_dict[key] = rescale_annData(
                    adata.copy(), list(adata_dict.items())[big_img][1]
                )
    for i, feat in enumerate(color):
        if same_scale:
            color_min = np.inf
            color_max = -np.inf
            for j, (title, adata) in enumerate(adata_dict.items()):
                if (
                    feat in adata.obs.columns
                    and adata.obs[feat].dtype.name != "category"
                ):
                    feat_min = np.min(adata.obs[feat])
                    feat_max = np.max(adata.obs[feat])
                elif feat in adata.var_names:
                    feat_min = np.min(adata[:, feat].X)
                    feat_max = np.max(adata[:, feat].X)
                else:
                    if feat not in adata.obs.columns:
                        print(f"Feature '{feat}' not found!")
                    feat_min = float("inf")
                    feat_max = float("-inf")
                color_min = min(color_min, feat_min)
                color_max = max(color_max, feat_max)

        for j, (tmp_key, adata) in enumerate(adata_dict.items()):
            ax = axes[i * ncols + j]
            ax.set_aspect("equal")
            # vmin = color_min if same_scale else None
            vmax = color_max if same_scale else None

            if "cmap" in kwargs.keys():
                if kwargs["cmap"] is not None and type(kwargs["cmap"]) == str:
                    kwargs["cmap"] = mpl.colormaps.get_cmap(kwargs["cmap"])
                    kwargs["cmap"].set_under(na_color)

            # Check if a variable in adata.obs is continuous or categorical
            if feat in adata.obs.columns and adata.obs[feat].dtype.name == "category":
                # Use the `sc.pl.embedding()` function to create the plot with the categorical variable
                ax.imshow(adata.uns["spatial"]["HE"])
                sc.pl.embedding(
                    adata, color=[feat], title=None, show=False, ax=ax, **kwargs
                )

            else:  # continuous variable
                # Create a divider for the existing axes instance
                divider = make_axes_locatable(ax)

                # Append axes to the right of ax
                cax = divider.append_axes(
                    "right", size=f"{colorbar_width}%", pad=colorbar_pad
                )

                sc.pl.embedding(
                    adata,
                    color=feat,
                    title=None,
                    show=False,
                    ax=ax,
                    vmax=vmax,
                    colorbar_loc=None,
                    **kwargs,
                )

                # manually add a colorbar
                fig.colorbar(ax.collections[0], cax=cax, orientation="vertical")

            # y-axis labels
            if j == 0:
                ax.set_ylabel(feat, fontsize="medium", fontstyle="italic")
            else:
                ax.set_ylabel(None)

            # x-axis labels
            ax.set_xlabel(None)
            if i == 0:
                ax.set_title(titles[j])
            else:
                ax.set_title(None)

    plt.tight_layout()

    if filename != None:
        plt.savefig(filename)

    plt.show()

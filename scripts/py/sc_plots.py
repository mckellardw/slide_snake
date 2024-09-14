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
    adata,
    xlim=[0, 20000],
    line_width=2,
    line_color="b",
    title="Knee plot",
    verbose=False,
):
    import matplotlib.pyplot as plt

    expected_num_cells = 10000

    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.loglog(range(len(knee)), knee, linewidth=line_width, color=line_color)
    #     ax.axvline(x=knee[expected_num_cells], linewidth=3, color="k")
    #     ax.axhline(y=expected_num_cells, linewidth=3, color="k")

    ax.set_ylabel("UMI Counts")
    ax.set_xlabel("Ranked Barcodes")
    ax.set_title(title)

    plt.xlim(xlim)
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


# scanpy version of seuListPlot - grids of plots for a single embedding
from mpl_toolkits.axes_grid1 import make_axes_locatable


# import mpl_toolkits.axes_grid1.inset_locator as mpl_il
def plot_grid_of_embeddings(
    adata_dict,
    color,
    ncols=None,
    titles=None,
    figsize=(10, 10),
    na_color="lightgrey",
    colorbar_width=10,
    colorbar_pad=0.05,
    same_scale=False,
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
                    kwargs["cmap"] = cm.get_cmap(kwargs["cmap"])
                    kwargs["cmap"].set_under(na_color)

            # Check if a variable in adata.obs is continuous or categorical
            if feat in adata.obs.columns and adata.obs[feat].dtype.name == "category":
                # Use the `sc.pl.embedding()` function to create the plot with the categorical variable
                sc.pl.embedding(
                    adata,
                    color=feat,
                    title=None,
                    show=False,
                    ax=ax,
                    legend_loc=None,
                    colorbar_loc=None,
                    **kwargs,
                )

                # Manually add a legend for the categorical variable
                categories = adata.obs[feat].cat.categories.tolist()
                legend_handles = [
                    plt.Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        markerfacecolor=color,
                        markersize=8,
                    )
                    for color in ax.collections[0].get_facecolor()
                ]

                # Create the legend
                ax.legend(
                    legend_handles,
                    categories,
                    title=feat,
                    fontsize="small",
                    bbox_to_anchor=(1, 0.0, (1 / colorbar_width), 1),
                    # loc='center left', bbox_to_anchor=(1, 0.5)
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

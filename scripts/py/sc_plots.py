# Plotting functions for use with scanpy
import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib import cm

# import anndata as ad


# Knee plot to quality check UMI counts for single-cell data
def knee_plot(
    adata,
    xlim=[0, 20000],
    line_width=2,
    line_color="b",
    title="Knee plot",
    expected_num_cells=10000,
    figsize=(10, 7),
    verbose=False,
):
    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]

    if expected_num_cells > len(knee):
        expected_num_cells = len(knee) - 1
    fig, ax = plt.subplots(figsize=figsize)

    ax.loglog(range(len(knee)), knee, linewidth=line_width, color=line_color)

    ax.axvline(x=expected_num_cells, linewidth=1, color="k")
    ax.axhline(y=knee[expected_num_cells], linewidth=1, color="k")

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



def ad_scatter(
    adata: ad.AnnData,
    x: str,
    y: str,
    x_layer: str = None,
    y_layer: str = None,
    color: str = None,
    alpha: float = 1,
    title: str = "Scatter Plot",
    figsize: tuple = (10, 8)
):
    """
    Create a scatter plot from an AnnData object. Compare 2 values.

    Parameters:
    -----------
    adata : anndata.AnnData
        The AnnData object containing the data.
    x : str
        The feature or metadata to plot on the x-axis.
    y : str
        The feature or metadata to plot on the y-axis.
    x_layer : str, optional
        The layer to use for the x-axis if x is a gene. Default is None (uses .X).
    y_layer : str, optional
        The layer to use for the y-axis if y is a gene. Default is None (uses .X).
    color : str, optional
        The feature or metadata to use for coloring the points.
    title : str, optional
        The title of the plot. Default is "Scatter Plot".
    figsize : tuple, optional
        The size of the figure. Default is (10, 8).

    Returns:
    --------
    fig : matplotlib.figure.Figure
        The created figure object.
    ax : matplotlib.axes.Axes
        The created axes object.
    """
    
    def get_feature_data(feature, layer):
        if feature in adata.var_names:
            if layer is None:
                return adata[:, feature].X.flatten()
            elif layer in adata.layers:
                return adata[:, feature].layers[layer].flatten()
            else:
                raise ValueError(f"Layer '{layer}' not found in the AnnData object.")
        elif feature in adata.obs.columns:
            return adata.obs[feature]
        else:
            raise ValueError(f"Feature '{feature}' not found in AnnData object.")

    x_data = get_feature_data(x, x_layer)
    y_data = get_feature_data(y, y_layer)

    fig, ax = plt.subplots(figsize=figsize)

    if color:
        color_data = get_feature_data(color, None)
        scatter = ax.scatter(
            x_data, 
            y_data, 
            c=color_data, 
            alpha=alpha,
            cmap='viridis'
        )
        plt.colorbar(
            scatter, 
            ax=ax, label=color)
    else:
        ax.scatter(
            x_data, 
            y_data,
            alpha=alpha
        )

    ax.set_xlabel(f"{x} ({'obs' if x in adata.obs.columns else f'var, layer: {x_layer if x_layer else "X"}'})")
    ax.set_ylabel(f"{y} ({'obs' if y in adata.obs.columns else f'var, layer: {y_layer if y_layer else "X"}'})")
    ax.set_title(title)

    plt.tight_layout()
    return fig, ax

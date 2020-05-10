import os
from typing import Tuple, Optional, Union
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from adjustText import adjust_text
from matplotlib.colorbar import ColorbarBase
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA
import functools
import warnings

from mspypeline.helpers import get_number_rows_cols_for_fig, plot_annotate_line, get_legend_elements, get_plot_name_suffix

FIG_FORMAT = ".pdf"


def linear(x, m, b):
    return m * x + b


def save_plot(name):
    """

    Parameters
    ----------
    name

    Returns
    -------

    """
    def decorator_save_plot(func):
        @functools.wraps(func)
        def wrapper_save_plot(*args, **kwargs):
            # get values with default
            plot_name = kwargs.get("name", name)
            save_path = kwargs.get("save_path", None)
            # suffix
            df_to_use = kwargs.get("df_to_use", None)
            level = kwargs.get("level", None)
            # file format
            fig_format = kwargs.get("fig_format", FIG_FORMAT)
            # run original function
            fig, ax = func(*args, **kwargs)
            # save if wanted
            if save_path is not None:
                try:
                    suffix = get_plot_name_suffix(df_to_use=df_to_use, level=level)
                    plot_name = f"{plot_name}{suffix}" + fig_format
                    plot_name = plot_name.format_map(kwargs)
                    res_path = os.path.join(save_path, plot_name)
                    fig.savefig(res_path, dpi=200, bbox_inches="tight")
                except PermissionError:
                    warnings.warn("Permission error in function %s. Did you forget to close the file?",
                                  str(func).split(" ")[1])
            return fig, ax
        return wrapper_save_plot
    return decorator_save_plot


_save_plot_decorator_doc = """
    save_path
    plot_name
    df_to_use
    level
    fig_format
"""


class QuantileNormalize(colors.Normalize):
    def __init__(self, quantiles: pd.Series):
        self.quantiles = quantiles
        assert self.quantiles.index.min() >= 0
        assert self.quantiles.index.max() <= 1
        self.outputs = np.concatenate((self.quantiles.index.values, [1]))
        super().__init__(0, 1)

    def __call__(self, value, clip=None):
        index = np.searchsorted(self.quantiles, value)
        return np.ma.array(self.outputs[index], fill_value=0)


def save_volcano_results(
        volcano_data: pd.DataFrame, unique_g1: pd.Series = None, unique_g2: pd.Series = None, g1: str = "group1",
        g2: str = "group2", col: str = "adjpval", intensity_label: str = "", save_path=".",
        show_suptitle: bool = True, fchange_threshold: float = 2, scatter_size: float = 10,
        n_labelled_proteins: int = 10
):
    """
    Saves multiple csv files and images containing the information of the volcano plot

    Parameters
    ----------
    volcano_data
        DataFrame containing data for the volcano plot with columns logFC and column specified under col. The index
        should be protein names or gene names
    unique_g1
        Series containing intensities of proteins or genes unique to group one
    unique_g2
        Series containing intensities of proteins or genes unique to group two
    g1
        Name of group one
    g2
        Name of group two
    col
        Column name containing p values
    intensity_label
        From which intensities were the fold changes calculated
    save_path
        path under which the results will be saved
    show_suptitle
        Should the figure title be shown
    fchange_threshold
        Minimum fold change threshold to be labelled significant
    scatter_size
        size of the points in the scatter plots
    n_labelled_proteins
        number of points that will be marked in th plot


    """
    plt.close("all")

    col_mapping = {"adjpval": "adjusted p value", "pval": "unadjusted p value"}

    if save_path is not None:
        # save all values
        volcano_data.to_csv(os.path.join(save_path,
                                      f"volcano_plot_data_{g1}_vs_{g2}_full_{col_mapping[col].replace(' ', '_')}.csv"))
        # save significant values
        volcano_data[volcano_data[col] < 0.05].to_csv(
            os.path.join(save_path,
                         f"volcano_plot_data_{g1}_vs_{g2}_significant_{col_mapping[col].replace(' ', '_')}.csv"))
        # save unique values
        if unique_g1 is not None:
            unique_g1.to_csv(os.path.join(save_path, f"volcano_plot_data_{g1}_vs_{g2}_unique_{g1}.csv"), header=True)
        if unique_g2 is not None:
            unique_g2.to_csv(os.path.join(save_path, f"volcano_plot_data_{g1}_vs_{g2}_unique_{g2}.csv"), header=True)

    def get_volcano_significances(fchange, pval, fchange_threshold):
        if pval > 0.05 or abs(fchange) < np.log2(fchange_threshold):
            return "ns"
        elif fchange >= 0:
            return "up"
        elif fchange < 0:
            return "down"
        else:
            raise ValueError(f"heisenbug: fold change: {fchange}, p value: {pval}")

    significance_to_color = {"ns": "gray", "up": "red", "down": "blue"}
    significance_to_label = {"ns": "non-significant", "up": f"upregulated in {g2}", "down": f"upregulated in {g1}"}

    # plot
    fig = plt.figure(figsize=(7, 7))

    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 8, 1])
    ax_unique_down: plt.Axes = plt.subplot(gs[0])
    ax: plt.Axes = plt.subplot(gs[1])
    ax_unique_up: plt.Axes = plt.subplot(gs[2])

    # hide the spines between ax and ax2
    ax_unique_down.spines['right'].set_visible(False)
    ax_unique_up.spines['left'].set_visible(False)
    ax_unique_down.yaxis.tick_left()
    ax_unique_up.yaxis.tick_right()
    ax_unique_up.yaxis.set_label_position("right")
    # hide the xticks
    ax_unique_down.tick_params(which='both', bottom=False, labelbottom=False)
    ax_unique_up.tick_params(which='both', bottom=False, labelbottom=False)

    # non sign gray, left side significant blue, right side red
    significances_plot = [get_volcano_significances(log_fold_change, p_val, fchange_threshold)
                          for log_fold_change, p_val in zip(volcano_data["logFC"], volcano_data[col])]
    for regulation in significance_to_color:
        mask = [x == regulation for x in significances_plot]
        ax.scatter(volcano_data["logFC"][mask], -np.log10(volcano_data[col])[mask], s=scatter_size,
                   color=significance_to_color[regulation], label=f"{sum(mask)} {significance_to_label[regulation]}")
    # get axis bounds for vertical and horizontal lines
    ymin, ymax = ax.get_ybound()
    xmin, xmax = ax.get_xbound()
    m = max(abs(xmin), xmax)
    # center the plot around 0
    ax.set_xlim(left=-1 * m, right=m)
    # update the x bounds
    xmin, xmax = ax.get_xbound()
    axline_kwargs = dict(linestyle="--", color="black", alpha=0.5, linewidth=1)
    # add line at significance threshold
    if any(volcano_data[col] < 0.05):
        x_offset = (np.log2(fchange_threshold) / xmax) / 2
        ax.axhline(-np.log10(0.05), **axline_kwargs, xmin=0, xmax=0.5 - x_offset)
        ax.axhline(-np.log10(0.05), **axline_kwargs, xmin=0.5 + x_offset, xmax=1)

    # add lines for minimum fold change threshold
    y_percentage = (-np.log10(0.05) + abs(ymin)) / (ymax + abs(ymin))
    if fchange_threshold > 0:
        ax.axvline(-np.log2(fchange_threshold), **axline_kwargs, ymin=y_percentage, ymax=1)
        ax.axvline(np.log2(fchange_threshold), **axline_kwargs, ymin=y_percentage, ymax=1)
    # plot unique values with mean intensity at over maximum
    ax_unique_down.scatter([0] * len(unique_g1), unique_g1, s=scatter_size, color="dodgerblue",
                           label=f"{len(unique_g1)} unique in {g1}")
    ax_unique_up.scatter([0] * len(unique_g2), unique_g2, s=scatter_size, color="coral",
                         label=f"{len(unique_g2)} unique in {g2}")
    # adjust bounds for unique axis
    ymin_down, ymax_down = ax_unique_down.get_ybound()
    ymin_up, ymax_up = ax_unique_up.get_ybound()
    ax_unique_down.set_ylim(bottom=min(ymin_down, ymin_up), top=max(ymax_down, ymax_up))
    ax_unique_up.set_ylim(bottom=min(ymin_down, ymin_up), top=max(ymax_down, ymax_up))

    # figure stuff
    if show_suptitle:
        fig.suptitle(f"{g1} vs {g2}")
    ax.set_xlabel(f"{intensity_label} Fold Change")
    ax.set_ylabel(r"-$Log_{10}$" + f" {col_mapping[col]}")
    ax_unique_down.set_ylabel(intensity_label)
    ax_unique_up.set_ylabel(intensity_label)
    fig.legend(bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False)
    res_path = os.path.join(save_path, f"volcano_{g1}_{g2}_no_annotation_{col_mapping[col].replace(' ', '_')}" + FIG_FORMAT)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(res_path, dpi=200, bbox_inches="tight")
    significant_upregulated = volcano_data[
        (volcano_data["logFC"] > np.log2(fchange_threshold)) & (volcano_data[col] < 0.05)
    ].sort_values(by=[col], ascending=True).head(n_labelled_proteins)
    significant_downregulated = volcano_data[
        (volcano_data["logFC"] < -np.log2(fchange_threshold)) & (volcano_data[col] < 0.05)
    ].sort_values(by=[col], ascending=True).head(n_labelled_proteins)
    significant = pd.concat([significant_upregulated, significant_downregulated])
    texts = []
    for log_fold_change, p_val, gene_name in zip(significant["logFC"], significant[col], significant.index):
        texts.append(ax.text(log_fold_change, -np.log10(p_val), gene_name, ha="center", va="center", fontsize=8))
    adjust_text(texts, arrowprops=dict(width=0.15, headwidth=0, color='gray', alpha=0.6), ax=ax)
    if save_path is not None:
        res_path = os.path.join(save_path, f"volcano_{g1}_{g2}_annotation_{col_mapping[col].replace(' ', '_')}" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")
    # TODO scatter plot of significant genes
    return fig, (ax, ax_unique_down, ax_unique_up)


@save_plot("pca_overview")
def save_pca_results(pca_data: pd.DataFrame, pca_fit: PCA = None, normalize: bool = True, save_path: Optional[str] = ".",
                     show_suptitle: bool = True, **kwargs):
    """
        Saves image containing the pca results

        Parameters
        ----------
        pca_data:
            DataFrame containing transformed/dimensionally-reduced data with which PCA was performed
        pca_fit:
            PCA object that was fitted to normalized input data
        normalize:
            Boolean whether the transformed data should be normalized with the singular values before plotting
        save_path:
            path under which the results will be saved
        show_suptitle:
            Should the figure title be shown


    """
    plt.close("all")
    n_components = pca_data.shape[0]
    singular_values = np.ones(n_components)
    color_map = {value: f"C{i}" for i, value in enumerate(pca_data.columns.get_level_values(0).unique())}
    color_map.update(kwargs.get("color_map", {}))
    if normalize and pca_fit is None:
        # TODO warning
        pass
    elif normalize and pca_fit is not None:
        singular_values = pca_fit.singular_values_
    fig, axarr = plt.subplots(n_components, n_components, figsize=(14, 14))
    for row in range(n_components):
        row_pc = row + 1
        for col in range(n_components):
            col_pc = col + 1
            if row > col:
                ax = axarr[col, row]
                ax.scatter(
                    pca_data.loc[f"PC_{row_pc}"] / singular_values[row],
                    pca_data.loc[f"PC_{col_pc}"] / singular_values[col],
                    c=[color_map.get(name, "blue") for name in pca_data.columns.get_level_values(0)])
                ax.set_xlabel(f"PC_{row_pc}")
                ax.set_ylabel(f"PC_{col_pc}")

    if show_suptitle:
        fig.suptitle(f"{kwargs['df_to_use']} intensity", fontsize="xx-large")
    legend_elements = get_legend_elements(labels=pca_data.columns.get_level_values(0).unique(), color_map=color_map)
    fig.legend(handles=legend_elements, bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False, fontsize=20)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, axarr


def save_pathway_analysis_results(
        protein_intensities: pd.DataFrame, significances: pd.DataFrame = None,
        pathway: str = "", show_suptitle: bool = False, level=0,
        threshold: float = 0.05, intensity_label: str = "", save_path=".", **kwargs
):
    plt.close("all")
    protein_intensities.to_csv(os.path.join(save_path, f"{pathway}_level_{level}_{intensity_label}_intensities.csv"))
    level_keys = list(protein_intensities.columns.get_level_values(0).unique())
    n_rows, n_cols = get_number_rows_cols_for_fig(protein_intensities.index)
    fig, axarr = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, int(n_rows * len(level_keys) / 1.5)))
    color_map = {value: f"C{i}" for i, value in enumerate(level_keys)}
    color_map.update(kwargs.get("color_map", {}))
    if show_suptitle:
        fig.suptitle(pathway)
    for protein, (pos, ax) in zip(protein_intensities.index, np.ndenumerate(axarr)):
        ax.scatter(protein_intensities.loc[protein], [level_keys.index(c) for c in protein_intensities.columns.get_level_values(0)],
                   c=[color_map[c] for c in protein_intensities.columns.get_level_values(0)])
        ax.set_title(protein)
        ax.set_ylim((-1, len(level_keys)))
        ax.set_yticks([i for i in range(len(level_keys))])
        ax.set_yticklabels(level_keys)
        ax.set_xlabel(intensity_label)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    if save_path is not None:
        res_path = os.path.join(save_path, f"{pathway}_level_{level}_{intensity_label}_no_labels" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")

    if significances is not None:
        significances.to_csv(os.path.join(save_path, f"{pathway}_level_{level}_{intensity_label}_pvalues.csv"))
        for protein, (pos, ax) in zip(protein_intensities.index, np.ndenumerate(axarr)):
            # adjust axis height based on number of significant differences
            to_annotate = significances.loc[protein]
            to_annotate = to_annotate[to_annotate <= threshold]
            xmin, xmax = ax.get_xbound()
            ax.set_xlim(right=xmax * (1 + to_annotate.shape[0] * 0.015))
            for i, (index, pval) in enumerate(to_annotate.items()):
                plot_annotate_line(ax, level_keys.index(index[0]), level_keys.index(index[1]), xmax * (1 + i * 0.015) - 0.005, pval)

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        if save_path is not None:
            res_path = os.path.join(save_path, f"{pathway}_level_{level}_{intensity_label}" + FIG_FORMAT)
            fig.savefig(res_path, dpi=200, bbox_inches="tight")
    return fig, axarr


@save_plot("boxplot")
def save_boxplot_results(
        protein_intensities: pd.DataFrame, intensity_label: str = "Intensity",
        plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, vertical: bool = False, **kwargs
):
    """
    Boxplot of intensities

    Parameters
    ----------
    protein_intensities
        DataFrame where each column are the intensities to boxplot, column names will be used as labels
    save_path
        path under which the results will be saved
    level
        level from with the data comes from. used for the save path
    intensity_label
        label of the x axis of the plot
    plot
        asd
    vertical

    kwargs
        accepts kwargs

    """
    # TODO give colors to the different groups
    if plot is None:
        plt.close("all")
        fig, ax = plt.subplots(figsize=(14, 1 + len(protein_intensities.columns) // 3))
    else:
        fig, ax = plot
    # indicate overall median with a line
    median_value = np.nanmedian(protein_intensities.values.flatten())
    line_kwargs = dict(color="black", alpha=0.5, linewidth=1)
    if vertical:
        ax.axhline(median_value, **line_kwargs)
    else:
        ax.axvline(median_value, **line_kwargs)
    # convert the data into a list of lists and filter nan values
    data = [
        protein_intensities.loc[~pd.isna(protein_intensities.loc[:, c]), c].tolist()
        for c in protein_intensities.columns
    ]
    ax.boxplot(data, vert=vertical, labels=protein_intensities.columns)
    if vertical:
        ax.set_ylabel(intensity_label)
    else:
        ax.set_xlabel(intensity_label)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


def save_relative_std_results(
        intensities: pd.DataFrame, name: str, df_to_use: str = "raw", intensity_label: str = "Intensity",
        bins=(10, 20, 30), save_path: Optional[str] = ".", cmap: dict = None, **kwargs
):
    """
    Relative standard deviations of passed intensities with color marking based on the specified bins and color map

    Parameters
    ----------
    intensities
        DataFrame with experiment intensities to be plotted
    name
        name of the overall experiment
    df_to_use
        from which intensities was the data gathered
    intensity_label
        name of the intensities for the x label
    bins
        in which bins should the standard deviations be categorized
    save_path
        path to which the figure will be saved
    cmap
        mapping for the digitized labels to a color
    kwargs
        accepts kwargs

    Returns
    -------
    figure and axis of the plot

    """
    # TODO add percentage to absolute numbers
    # TODO see if code could be optimized
    plt.close("all")

    bins = np.array(bins)
    if "log2" in df_to_use:
        bins = np.log2(bins)

    cm = {0: "navy", 1: "royalblue", 2: "skyblue", 3: "darkgray"}
    if cmap is not None:
        cm.update(cmap)

    relative_std_percent = intensities.std(axis=1) / intensities.mean(axis=1) * 100

    inds = np.digitize(relative_std_percent, bins).astype(int)
    plot_colors = pd.Series([cm.get(x, "black") for x in inds], index=relative_std_percent.index)
    color_counts = {color: (plot_colors == color).sum() for color in plot_colors.unique()}

    fig, ax = plt.subplots(1, 1, figsize=(14, 7))
    ax.scatter(intensities.mean(axis=1), relative_std_percent, c=plot_colors, marker="o", s=(2 * 72. / fig.dpi) ** 2,
               alpha=0.8)
    ax.set_xlabel(f"Mean {intensity_label}")
    ax.set_ylabel("Relative Standard deviation [%]")
    if "log2" not in df_to_use:
        ax.set_xscale('log')
    xmin, xmax = ax.get_xbound()
    cumulative_count = 0
    for i, bin_ in enumerate(bins):
        cumulative_count += color_counts.get(cm[i], 0)
        ax.axhline(bin_, color=cm[i])
        ax.text(xmin, bin_, cumulative_count)

    if save_path is not None:
        res_path = os.path.join(save_path, f"rel_std_{name}_{df_to_use}" + FIG_FORMAT)
        fig.savefig(res_path, dpi=200, bbox_inches="tight")
        plt.close(fig)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("detected_counts")
def save_detection_counts_results(
        counts: pd.DataFrame, intensity_label: str = "Intensity", show_suptitle: bool = True, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    f"""

    Parameters
    ----------
    counts
        DataFrame containing the counts to be plotted
    intensity_label
        label of the dataframe
    show_suptitle
        should the figure title be shown
    {_save_plot_decorator_doc}
    kwargs
        accepts kwargs

    Returns
    -------
    figure and axis of the plot

    """
    plt.close("all")

    n_rows_experiment, n_cols_experiment = get_number_rows_cols_for_fig(counts.columns)
    fig, axarr = plt.subplots(n_rows_experiment, n_cols_experiment, squeeze=True,
                              figsize=(5 * n_cols_experiment, 3 * n_rows_experiment))
    if show_suptitle:
        fig.suptitle(f"Detection counts from {intensity_label}")

    for (pos, ax), col in zip(np.ndenumerate(axarr), counts.columns):
        col_data = counts.loc[:, col]
        col_data = col_data[~pd.isna(col_data)]

        ax.set_title(f"{col},\ntotal detected: {int(col_data.sum())}")
        ax.barh(col_data.index, col_data, color="skyblue")
        for y, value in zip(col_data.index, col_data):
            ax.text(col_data.max() / 2, y, value,
                    verticalalignment='center', horizontalalignment='center')

        ax.set_yticks(col_data.index)
        ax.set_yticklabels([f"detected in {i} replicates" for i in col_data.index])
        ax.set_xlabel("Counts")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, axarr


@save_plot("kde")
def save_kde_results(
        intensities: pd.DataFrame, quantile_range: Optional[np.array] = None, n_points: int = 1000,
        cmap: Union[str, colors.Colormap] = "viridis", plot: Optional[Tuple[plt.Figure, plt.Axes]] = None,
        intensity_label: str = "Intensity", **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    if plot is not None:
        fig, ax = plot
    else:
        plt.close("all")
        fig, ax = plt.subplots(1, 1)

    if quantile_range is None:
        quantile_range = np.arange(0.05, 1, 0.05)

    for col in intensities.columns:
        intensity_quantiles = intensities.loc[:, col].quantile(quantile_range)
        kde_fit = gaussian_kde(intensities.loc[~pd.isna(intensities.loc[:, col]), col])

        x = np.linspace(intensities.loc[:, col].min() * 0.9, intensities.loc[:, col].max() * 1.1, n_points)
        y = kde_fit.evaluate(x)
        # Create a set of line segments so that we can color them individually
        # This creates the points as a N x 1 x 2 array so that we can stack points
        # together easily to get the segments. The segments array for line collection
        # needs to be (numlines) x (points per line) x 2 (for x and y)
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        #
        norm = QuantileNormalize(quantiles=intensity_quantiles)
        lc = LineCollection(segments, cmap=cmap, norm=norm, alpha=0.5)
        # Set the values used for colormapping
        lc.set_array(x)
        lc.set_linewidth(1)
        line = ax.add_collection(lc)

    ax.set_ylabel("Density")
    ax.set_xlabel(intensity_label)

    ax.autoscale_view()

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("n_proteins_vs_quantile")
def save_n_proteins_vs_quantile_results(
        quantiles: pd.DataFrame, n_proteins: pd.Series, nstd: int = 1, cmap: Union[str, colors.Colormap] = "viridis",
        plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, cbar_ax: Optional[plt.Axes] = None, intensity_label: str = "Intensity",
        fill_between: bool = False, **kwargs
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    if plot is not None:
        fig, ax = plot
    else:
        plt.close("all")
        fig, ax = plt.subplots(1, 1, figsize=(14, 7))

    if not isinstance(cmap, colors.Colormap):
        cmap = cm.get_cmap(cmap)

    m = n_proteins.sort_values()
    for quant in quantiles.index:
        ax.scatter(quantiles.loc[quant, :], n_proteins, c=[cmap(quant)] * len(n_proteins), alpha=0.5)
        popt, pcov = curve_fit(linear, n_proteins, quantiles.loc[quant, :])
        fit = linear(m, *popt)

        ax.plot(fit, m, color=cmap(quant))

        if fill_between:
            perr = np.sqrt(np.diag(pcov))
            popt_up = popt + nstd * perr
            popt_dw = popt - nstd * perr
            fit_up = linear(m, *popt_up)
            fit_dw = linear(m, *popt_dw)
            ax.fill_betweenx(m, fit_dw, fit_up, alpha=.15, color=cmap(quant))

    if cbar_ax is None:
        cbar_ax = fig.add_axes([0.2, 0.00, 0.6, 0.02])  # [left, bottom, width, height]
    cb = ColorbarBase(cbar_ax, cmap=cmap, orientation="horizontal")
    cb.set_label("quantile")

    ax.set_ylabel("# detected proteins")
    ax.set_xlabel(intensity_label)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, (ax, cbar_ax)


@save_plot("normalization_overview")
def save_normalization_overview_results(
        quantiles, n_proteins, intensities, protein_intensities,
        height: int = 15, intensity_label: str = "Intensity", **kwargs
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes, plt.Axes, plt.Axes]]:
    f"""
    
    Parameters
    ----------
    quantiles
    n_proteins
    intensities
    protein_intensities
    height
    intensity_label
    {_save_plot_decorator_doc}
    kwargs

    Returns
    -------

    """
    fig = plt.figure(figsize=(18, 18))
    gs = fig.add_gridspec(height, 2)
    ax_density = fig.add_subplot(gs[0:height // 2, 0])
    ax_nprot = fig.add_subplot(gs[height // 2:height - 1, 0])
    ax_colorbar = fig.add_subplot(gs[height - 1, 0])
    ax_boxplot = fig.add_subplot(gs[0:height, 1])

    # order the boxplot data after the number of identified peptides
    boxplot_data = protein_intensities.loc[:, n_proteins.sort_values(ascending=False).index[::-1]]

    plot_kwargs = dict(intensity_label=intensity_label)
    plot_kwargs.update(**kwargs)
    plot_kwargs["save_path"] = None
    save_kde_results(intensities=intensities,  plot=(fig, ax_density), **plot_kwargs)
    save_n_proteins_vs_quantile_results(quantiles=quantiles, n_proteins=n_proteins, plot=(fig, ax_nprot), cbar_ax=ax_colorbar, **plot_kwargs)
    save_boxplot_results(boxplot_data, plot=(fig, ax_boxplot), vertical=False, **plot_kwargs)
    ax_density.set_xlim(ax_nprot.get_xlim())

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, (ax_nprot, ax_density, ax_colorbar, ax_boxplot)


@save_plot("intensities_heatmap")
def save_intensities_heatmap_result(
        intensities: pd.DataFrame, cmap: Union[str, colors.Colormap] = "autumn_r", cmap_bad="dimgray",
        cax: plt.Axes = None, plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, **kwargs
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    f"""
    
    Parameters
    ----------
    intensities
    cmap
    cmap_bad
    cax
    plot
    {_save_plot_decorator_doc}
    kwargs

    Returns
    -------

    """
    if plot is not None:
        fig, ax = plot
    else:
        plt.close("all")
        fig, ax = plt.subplots(figsize=(14, 16))  # TODO scale with intensities.shape

    if not isinstance(cmap, colors.Colormap):
        cmap = cm.get_cmap(cmap)
    if cmap_bad is not None:
        cmap.set_bad(color='dimgray')

    im = ax.imshow(intensities.values.T, aspect="auto", cmap=cmap)
    if cax is None:
        cbar = ax.figure.colorbar(im, ax=ax)
    else:
        cbar = ax.figure.colorbar(im, cax=cax)

    y_lim = ax.get_ylim()
    ax.set_yticks(np.linspace(0, len(intensities.columns) - 1, len(intensities.columns)))
    ax.set_yticklabels(intensities.columns)
    ax.set_ylim(*y_lim)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, (ax, cbar.ax)

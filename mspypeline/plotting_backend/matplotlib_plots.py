import os
from itertools import combinations
from typing import Tuple, Optional, Union, Callable, Dict, Iterable
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection
from adjustText import adjust_text
from matplotlib.colorbar import ColorbarBase
from matplotlib_venn import venn2, venn3
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
from scipy import stats
from sklearn.decomposition import PCA
import functools
import warnings

from mspypeline.helpers import get_number_rows_cols_for_fig, plot_annotate_line, get_legend_elements, \
    get_plot_name_suffix, get_intersection_and_unique, venn_names

FIG_FORMAT = ".pdf"


def linear(x, m, b):
    return m * x + b


def collect_plots_to_pdf(path: str, *args, dpi: int = 200):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if not path.endswith(".pdf"):
        path += ".pdf"
    with PdfPages(path) as pdf:
        for plot in args:
            figure = None
            if callable(plot):
                figure, axes = plot()
            elif isinstance(plot, Iterable):
                for x in plot:
                    if isinstance(x, plt.Figure):
                        figure = x
                        break
            elif isinstance(plot, plt.Figure):
                figure = plot
            if figure is not None:
                pdf.savefig(figure, dpi=dpi)
                plt.close(figure)


_get_path_and_name_kwargs_doc = """
Uses following kwargs to generate a path and file name.
save_path
df_to_use
level

Additionally, all parts of the passed name enclosed by "{}"
will be replaced from passed kwargs.
E.g. if the passed name is "plot_{stuff}", a stuff kwarg is expected and will be replaced by the passed value.
"""


def get_path_and_name_from_kwargs(name: str, **kwargs) -> Tuple[str, str]:
    f"""
    generates stuff

    Parameters
    ----------
    name
    kwargs
        {_get_path_and_name_kwargs_doc}

    Returns
    -------
    The path where the file will be saved and a name for the file without a file type.
    """
    save_path = kwargs.get("save_path", None)
    # get values with default
    # suffix
    df_to_use = kwargs.get("df_to_use", None)
    level = kwargs.get("level", None)

    suffix = get_plot_name_suffix(df_to_use=df_to_use, level=level)
    name = f"{name}{suffix}"
    name = name.format_map(kwargs)
    return save_path, name


def save_plot_func(
        fig: plt.Figure, path: str, plot_name: str, func: Callable, fig_format: str = FIG_FORMAT,
        dpi: int = 200, **kwargs
) -> None:
    """

    Parameters
    ----------
    fig
    path
    plot_name
    func
    fig_format
    dpi
    kwargs

    Returns
    -------

    """
    if path is not None:
        try:
            os.makedirs(path, exist_ok=True)
            res_path = os.path.join(path, plot_name)
            fig.savefig(res_path + fig_format, dpi=dpi, bbox_inches="tight")
        except PermissionError:
            warnings.warn(f"Permission error in function {str(func).split(' ')[1]}. Did you forget to close the file?")


def save_plot(plot_name: str):
    """
    Decorator to save figures, which are returned by the decorated function. Assumes that a tuple of figure, axes is
    returned. Plot is saved by using get_path_and_name_from_kwargs and save_plot_func.

    Parameters
    ----------
    plot_name

    Returns
    -------
    figure and axes of the decorated function

    """
    def decorator_save_plot(func):
        @functools.wraps(func)
        def wrapper_save_plot(*args, **kwargs):
            # run original function
            ret = func(*args, **kwargs)
            if ret is not None:
                path, pn = get_path_and_name_from_kwargs(name=plot_name, **kwargs)
                save_plot_func(ret[0], path, pn, func, **kwargs)
            return ret
        return wrapper_save_plot
    return decorator_save_plot


def save_csvs(name_map: Dict[str, str]):
    """
    Saves all dataframes as csv which are specified as dict keys. The values are the file names.

    Parameters
    ----------
    name_map

    Returns
    -------
    original function
    """
    def decorator_save_csvs(func):
        @functools.wraps(func)
        def wrapper_save_csvs(*args, **kwargs):
            for kwarg_name, file_name in name_map.items():
                df = kwargs.get(kwarg_name, None)
                if df is not None:
                    save_path, csv_name = get_path_and_name_from_kwargs(file_name, **kwargs)
                    os.makedirs(save_path, exist_ok=True)
                    if save_path is not None:
                        os.makedirs(os.path.dirname(save_path), exist_ok=True)
                        df.to_csv(os.path.join(save_path, csv_name) + ".csv", header=True)
            return func(*args, **kwargs)
        return wrapper_save_csvs
    return decorator_save_csvs


def save_venn_to_txt(name_map: Dict[str, str]):
    def decorator_save_venn(func):
        @functools.wraps(func)
        def wrapper_save_venn(*args, **kwargs):
            for kwarg_name, file_name in name_map.items():
                named_sets = kwargs.get(kwarg_name, None)
                if named_sets is not None:
                    if len(named_sets) > 6:
                        warnings.warn(f"Skipping save_venn_to_txt because more than 6 experiments were passed at once. "
                                      f"({len(named_sets)}")
                        continue
                    save_path, txt_name = get_path_and_name_from_kwargs(file_name, **kwargs)
                    for intersected, unioned, result in venn_names(named_sets):
                        # create name based on the intersections and unions that were done
                        intersected_name = "&".join(sorted(intersected))
                        unioned_name = "-" + "-".join(sorted(unioned)) if unioned else ""
                        res_path = os.path.join(save_path, f"{txt_name}_{intersected_name}{unioned_name}.txt")
                        # write all names line by line into the file
                        with open(res_path, "w") as out:
                            for re in result:
                                out.write(re + "\n")
            return func(*args, **kwargs)
        return wrapper_save_venn
    return decorator_save_venn


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


@save_csvs({"unique_g1": "volcano_plot_data_{g1}_vs_{g2}_unique_{g1}",
            "unique_g2": "volcano_plot_data_{g1}_vs_{g2}_unique_{g2}"})
def save_volcano_results(
        volcano_data: pd.DataFrame, unique_g1: pd.Series = None, unique_g2: pd.Series = None, g1: str = "group1",
        g2: str = "group2", col: str = "adjpval", intensity_label: str = "Intensity",
        show_suptitle: bool = True, fchange_threshold: float = 2, scatter_size: float = 10,
        n_labelled_proteins: int = 10, **kwargs
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes, plt.Axes]]:
    f"""
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
    show_suptitle
        Should the figure title be shown
    fchange_threshold
        Minimum fold change threshold to be labelled significant
    scatter_size
        size of the points in the scatter plots
    n_labelled_proteins
        number of points that will be marked in th plot
    kwargs
        {_get_path_and_name_kwargs_doc}


    """
    plt.close("all")

    col_mapping = {"adjpval": "adjusted p value", "pval": "unadjusted p value"}

    save_path, csv_name = get_path_and_name_from_kwargs("volcano_plot_data_{g1}_vs_{g2}_full_{p}",
                                                        g1=g1, g2=g2, p=col_mapping[col].replace(' ', '_'), **kwargs)
    if save_path is not None:
        volcano_data.to_csv(os.path.join(save_path, csv_name) + ".csv", header=True)
    save_path, csv_name = get_path_and_name_from_kwargs("volcano_plot_data_{g1}_vs_{g2}_significant_{p}",
                                                        g1=g1, g2=g2, p=col_mapping[col].replace(' ', '_'), **kwargs)
    if save_path is not None:
        volcano_data[volcano_data[col] < 0.05].to_csv(os.path.join(save_path, csv_name) + ".csv", header=True)

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
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    # save intermediate results
    path, plot_name = get_path_and_name_from_kwargs(name="volcano_{g1}_{g2}_no_annotation_{p}", g1=g1, g2=g2,
                                                         p=col_mapping[col].replace(' ', '_'), **kwargs)
    save_plot_func(fig, path, plot_name, save_volcano_results, **kwargs)

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

    path, plot_name = get_path_and_name_from_kwargs(name="volcano_{g1}_{g2}_annotation_{p}", g1=g1, g2=g2,
                                                         p=col_mapping[col].replace(' ', '_'), **kwargs)
    save_plot_func(fig, path, plot_name, save_volcano_results, **kwargs)
    # TODO scatter plot of significant genes
    return fig, (ax, ax_unique_down, ax_unique_up)


@save_plot("pca_overview")
def save_pca_results(
        pca_data: pd.DataFrame, pca_fit: PCA = None, normalize: bool = True, intensity_label: str = "Intensity",
        color_map: Optional[dict] = None, show_suptitle: bool = True, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
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
    intensity_label

    color_map

    show_suptitle:
        Should the figure title be shown


    """
    plt.close("all")
    n_components = pca_data.shape[0]
    singular_values = np.ones(n_components)
    base_color_map = {value: f"C{i}" for i, value in enumerate(pca_data.columns.get_level_values(0).unique())}
    color_map = {} if color_map is None else color_map
    base_color_map.update(color_map)
    if normalize and pca_fit is None:
        warnings.warn("Normalizing not possible when pca_fit is None")
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
                    c=[base_color_map.get(name, "blue") for name in pca_data.columns.get_level_values(0)])
                ax.set_xlabel(f"PC_{row_pc}")
                ax.set_ylabel(f"PC_{col_pc}")

    if show_suptitle:
        fig.suptitle(intensity_label, fontsize="xx-large")
    legend_elements = get_legend_elements(labels=pca_data.columns.get_level_values(0).unique(), color_map=base_color_map)
    fig.legend(handles=legend_elements, bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False, fontsize=20)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, axarr


@save_csvs({"protein_intensities": "{pathway}_protein_intensities",
            "significances": "{pathway}_pvalues"})
def save_pathway_analysis_results(
        protein_intensities: pd.DataFrame, significances: pd.DataFrame = None, pathway: str = "",
        show_suptitle: bool = False, threshold: float = 0.05, intensity_label: str = "Intensity", **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    f"""
    
    Parameters
    ----------
    protein_intensities
    significances
    pathway
    show_suptitle
    threshold
    intensity_label
    kwargs
        {_get_path_and_name_kwargs_doc}

    Returns
    -------

    """
    plt.close("all")
    level_keys = list(protein_intensities.columns.get_level_values(0).unique())
    n_rows, n_cols = get_number_rows_cols_for_fig(protein_intensities.index)
    fig, axarr = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, int(n_rows * len(level_keys) / 1.5)))
    color_map = {value: f"C{i}" for i, value in enumerate(level_keys)}
    color_map.update(kwargs.get("color_map", {}))  # TODO move to function deceleration
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

    path, plot_name = get_path_and_name_from_kwargs(name="{pathway}_no_labels", pathway=pathway, **kwargs)
    save_plot_func(fig, path, plot_name, save_pathway_analysis_results, **kwargs)

    if significances is not None:
        for protein, (pos, ax) in zip(protein_intensities.index, np.ndenumerate(axarr)):
            # adjust axis height based on number of significant differences
            to_annotate = significances.loc[protein]
            to_annotate = to_annotate[to_annotate <= threshold]
            xmin, xmax = ax.get_xbound()
            ax.set_xlim(right=xmax * (1 + to_annotate.shape[0] * 0.015))
            for i, (index, pval) in enumerate(to_annotate.items()):
                plot_annotate_line(ax, level_keys.index(index[0]), level_keys.index(index[1]), xmax * (1 + i * 0.015) - 0.005, pval)

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        path, plot_name = get_path_and_name_from_kwargs(name="{pathway}", pathway=pathway, **kwargs)
        save_plot_func(fig, path, plot_name, save_pathway_analysis_results, **kwargs)
    return fig, axarr


@save_plot("boxplot")
def save_boxplot_results(
        protein_intensities: pd.DataFrame, intensity_label: str = "Intensity",
        plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, vertical: bool = False, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    f"""
    Boxplot of intensities

    Parameters
    ----------
    protein_intensities
        DataFrame where each column are the intensities to boxplot, column names will be used as labels
    intensity_label
        label of the x axis of the plot
    plot
        asd
    vertical

    kwargs
        {_get_path_and_name_kwargs_doc}

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


@save_plot("rel_std_{experiment_name}")
def save_relative_std_results(
        intensities: pd.DataFrame, experiment_name: str, intensity_label: str = "Intensity",
        show_suptitle: bool = False, bins=(10, 20, 30), cmap: dict = None, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    f"""
    Relative standard deviations of passed intensities with color marking based on the specified bins and color map

    Parameters
    ----------
    intensities
        DataFrame with experiment intensities to be plotted
    experiment_name
        name of the overall experiment
    intensity_label
        name of the intensities for the x label
    show_suptitle
        should figure suptitle be shown
    bins
        in which bins should the standard deviations be categorized
    cmap
        mapping for the digitized labels to a color
    kwargs
        {_get_path_and_name_kwargs_doc}

    Returns
    -------
    figure and axis of the plot

    """
    # TODO add percentage to absolute numbers
    # TODO see if code could be optimized
    plt.close("all")

    bins = np.array(bins)
    if "Log_2" in intensity_label:
        bins = np.log2(bins)

    default_cm = {0: "navy", 1: "royalblue", 2: "skyblue", 3: "darkgray"}
    if cmap is not None:
        default_cm.update(cmap)

    relative_std_percent = intensities.std(axis=1) / intensities.mean(axis=1) * 100

    inds = np.digitize(relative_std_percent, bins).astype(int)
    plot_colors = pd.Series([default_cm.get(x, "black") for x in inds], index=relative_std_percent.index)
    color_counts = {color: (plot_colors == color).sum() for color in plot_colors.unique()}

    fig, ax = plt.subplots(1, 1, figsize=(14, 7))
    ax.scatter(intensities.mean(axis=1), relative_std_percent, c=plot_colors, marker="o", s=(2 * 72. / fig.dpi) ** 2,
               alpha=0.8)
    if show_suptitle:
        fig.suptitle(experiment_name)
    ax.set_xlabel(f"Mean {intensity_label}")
    ax.set_ylabel("Relative Standard deviation [%]")
    if "Log_2" not in intensity_label:
        ax.set_xscale('log')
    xmin, xmax = ax.get_xbound()
    cumulative_count = 0
    for i, bin_ in enumerate(bins):
        cumulative_count += color_counts.get(default_cm[i], 0)
        ax.axhline(bin_, color=default_cm[i])
        ax.text(xmin, bin_, cumulative_count)

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
    kwargs
        {_get_path_and_name_kwargs_doc}

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
    f"""
    
    Parameters
    ----------
    intensities
    quantile_range
    n_points
    cmap
    plot
    intensity_label
    kwargs
        {_get_path_and_name_kwargs_doc}

    Returns
    -------

    """
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
        plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, cbar_ax: Optional[plt.Axes] = None,
        intensity_label: str = "Intensity", fill_between: bool = False, **kwargs
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    f"""
    
    Parameters
    ----------
    quantiles
    n_proteins
    nstd
    cmap
    plot
    cbar_ax
    intensity_label
    fill_between
    kwargs
        {_get_path_and_name_kwargs_doc}

    Returns
    -------

    """
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
    kwargs
        {_get_path_and_name_kwargs_doc}

    Returns
    -------

    """
    fig = plt.figure(figsize=(18, 18), constrained_layout=True)
    gs = fig.add_gridspec(height, 2)
    ax_density = fig.add_subplot(gs[0:height // 2, 0])
    ax_nprot = fig.add_subplot(gs[height // 2:height - 1, 0])
    ax_colorbar = fig.add_subplot(gs[height - 1, 0])
    ax_boxplot = fig.add_subplot(gs[0:height, 1])

    fig.suptitle(f"Normalization overview for {intensity_label}")
    # order the boxplot data after the number of identified peptides
    boxplot_data = protein_intensities.loc[:, n_proteins.sort_values(ascending=False).index[::-1]]

    plot_kwargs = dict(intensity_label=intensity_label)
    plot_kwargs.update(**kwargs)
    plot_kwargs["save_path"] = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # ignore warning for mixing constrained layout with thight_layout
        save_kde_results(intensities=intensities,  plot=(fig, ax_density), **plot_kwargs)
        save_n_proteins_vs_quantile_results(quantiles=quantiles, n_proteins=n_proteins, plot=(fig, ax_nprot), cbar_ax=ax_colorbar, **plot_kwargs)
        save_boxplot_results(boxplot_data, plot=(fig, ax_boxplot), vertical=False, **plot_kwargs)
    ax_density.set_xlim(ax_nprot.get_xlim())

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, (ax_nprot, ax_density, ax_colorbar, ax_boxplot)


@save_plot("intensities_heatmap")
def save_intensities_heatmap_result(
        intensities: pd.DataFrame, cmap: Union[str, colors.Colormap] = "autumn_r", cmap_bad="dimgray",
        cax: plt.Axes = None, plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, vmax: Optional[float] = None,
        vmin: Optional[float] = None,
        intensity_label: str = "Intensity", show_suptitle: bool = True, **kwargs
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    f"""
    
    Parameters
    ----------
    intensities
    cmap
    cmap_bad
    cax
    plot
    vmax
    vmin
    intensity_label
    show_suptitle
    kwargs
        {_get_path_and_name_kwargs_doc}

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

    im = ax.imshow(intensities.values.T, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    if cax is None:
        cbar = ax.figure.colorbar(im, ax=ax)
    else:
        cbar = ax.figure.colorbar(im, cax=cax)

    if show_suptitle:
        fig.suptitle(f"Present or absent in {intensity_label}")
    ax.set_xlabel("Proteins")

    y_lim = ax.get_ylim()
    ax.set_yticks(np.linspace(0, len(intensities.columns) - 1, len(intensities.columns)))
    ax.set_yticklabels(intensities.columns)
    ax.set_ylim(*y_lim)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, (ax, cbar.ax)


@save_plot("detection_per_replicate")
def save_number_of_detected_proteins_results(
        all_heights: Dict[str, pd.Series], intensity_label: str = "Intensity", show_suptitle: bool = True, **kwargs
):
    plt.close("all")
    # determine number of rows and columns in the plot based on the number of experiments
    n_rows_experiment, n_cols_experiment = get_number_rows_cols_for_fig(all_heights.keys())
    fig, axarr = plt.subplots(n_rows_experiment, n_cols_experiment,
                              figsize=(5 * n_cols_experiment, 3 * n_rows_experiment))
    if show_suptitle:
        fig.suptitle(f"Number of detected proteins from {intensity_label}")

    for experiment, (pos, ax) in zip(all_heights.keys(), np.ndenumerate(axarr)):
        experiment_heights = all_heights[experiment]
        mean_height = experiment_heights[1:].mean()
        y_pos = [x for x in range(len(experiment_heights))]
        ax.barh(y_pos, experiment_heights, color="skyblue")

        for y, value in zip(y_pos, experiment_heights):
            ax.text(experiment_heights[0] / 2, y, value,
                    verticalalignment='center', horizontalalignment='center')
        ax.set_title(experiment)
        ax.axvline(mean_height, linestyle="--", color="black", alpha=0.6)
        ax.set_yticks([i for i in range(len(experiment_heights.index))])
        ax.set_yticklabels(experiment_heights.index)
        ax.set_xlabel("Counts")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, axarr


@save_plot("intensity_histograms")
def save_intensity_histogram_results(
        hist_data: pd.DataFrame, intensity_label: str = "Intensity", show_suptitle: bool = False,
        compare_to_remaining: bool = False, n_bins: int = 25, histtype="bar", color=None,
        plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, **kwargs
):
    if plot is not None:
        fig, axarr = plot
    else:
        plt.close("all")
        n_rows, n_cols = get_number_rows_cols_for_fig(hist_data.columns.get_level_values(0).unique())
        fig, axarr = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))

    if show_suptitle:
        fig.suptitle(f"{intensity_label} histograms")

    for col, (pos, ax) in zip(hist_data.columns.get_level_values(0).unique(), np.ndenumerate(axarr)):
        intensities = hist_data[col]
        try:
            labels = intensities.columns
        except AttributeError:
            labels = col

        if "Log_2" in intensity_label:
            bins = np.linspace(np.nanmin(intensities.values), np.nanmax(intensities.values), n_bins)
        else:
            bins = np.logspace(np.log2(np.nanmin(intensities.values)), np.log2(np.nanmax(intensities.values)), n_bins, base=2)

        ax.set_title(col)
        ax.hist(intensities.T, bins=bins, histtype=histtype, label=labels, color=color)

        if compare_to_remaining:
            remaining = hist_data.drop(col, axis=1)
            remaining = remaining.mean(axis=1)
            # remaining = remaining[intensities.notna()]
            ax.hist(remaining, bins=bins, histtype="step", alpha=0.5)
        if "Log_2" not in intensity_label:
            ax.set_xscale("log", basex=2)
        ax.set_xlabel(intensity_label)
        ax.set_ylabel("Counts")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, axarr


@save_plot("{full_name}_scatter")
def save_scatter_replicates_results(
        scatter_data: pd.DataFrame, intensity_label: str = "Intensity", show_suptitle: bool = False, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    plt.close("all")
    fig, ax = plt.subplots(1, 1, figsize=(7, 7))

    if show_suptitle:
        fig.suptitle()

    for rep1, rep2 in combinations(scatter_data.columns, 2):
        x1 = scatter_data.loc[:, rep1]
        x2 = scatter_data.loc[:, rep2]
        corr_mask = np.logical_and(x1.notna(), x2.notna())
        plot_mask = np.logical_or(x1.notna(), x2.notna())
        exp = r"$r^{2}$"
        ax.scatter(x1.fillna(x2.min() * 0.95)[plot_mask], x2.fillna(x2.min() * 0.95)[plot_mask],
                   label=f"{rep1} vs {rep2}, "
                         fr"{exp}: {stats.pearsonr(x1[corr_mask], x2[corr_mask])[0] ** 2:.4f}",
                   alpha=0.5, marker=".")
        ax.set_xlabel(intensity_label)
        ax.set_ylabel(intensity_label)

    fig.legend(frameon=False)
    if "Log_2" not in intensity_label:
        ax.set_xscale("log")
        ax.set_yscale("log")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("{full_name}_rank")
def save_rank_results(
        rank_data: pd.Series, interesting_proteins, intensity_label: str = "Intensity", full_name="Experiment",
        show_suptitle: bool = False, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    plt.close("all")
    if interesting_proteins.values():
        all_pathway_proteins = set.union(*(set(x) for x in interesting_proteins.values()))
    else:
        all_pathway_proteins = set()
    # get the experiment intensities calculate mean intensity for the experiment and sort from highest to lowest
    # TODO apply filter for rare proteins before here?
    # protein ranks vs intensity
    # create dict to map each protein its respective rank and mean intensity
    dic = {idx: (i, value) for i, (idx, value) in enumerate(rank_data.items())}

    found_proteins = set(rank_data.index)
    # get all proteins that are not part of any pathway
    non_pathway_proteins = found_proteins - all_pathway_proteins
    # get all proteins that are part of any pathway
    pathway_proteins = found_proteins & all_pathway_proteins
    rank_identified_proteins = [dic[protein][0] for protein in pathway_proteins]
    # plot the non pathway proteins
    x = [dic[protein][0] for protein in non_pathway_proteins]
    y = [dic[protein][1] for protein in non_pathway_proteins]

    fig, ax = plt.subplots(1, 1, figsize=(14, 7))
    ax.scatter(x, y, c=f"darkgray", s=10, alpha=0.3, marker=".", label="no pathway")
    # plot all proteins of a specific pathway
    for i, (pathway, proteins) in enumerate(interesting_proteins.items()):
        proteins = set(proteins) & found_proteins
        x = [dic[protein][0] for protein in proteins]
        y = [dic[protein][1] for protein in proteins]
        ax.scatter(x, y, c=f"C{i}", s=80, alpha=0.6, marker=".", label=pathway)

    # only if more than 0 proteins are identified
    if rank_identified_proteins:
        median_pathway_rank = int(np.median(rank_identified_proteins))
        median_intensity = rank_data.iloc[median_pathway_rank]
        xmin, xmax = ax.get_xbound()
        xm = (median_pathway_rank + abs(xmin)) / (abs(xmax) + abs(xmin))
        ymin, ymax = ax.get_ybound()
        ym = (median_intensity - ymin) / (ymax - ymin)
        # plot the median rank and intensity at that rank
        ax.axvline(median_pathway_rank, ymax=ym, linestyle="--", color="black", alpha=0.6)
        ax.axhline(median_intensity, xmax=xm, linestyle="--", color="black", alpha=0.6)
        ax.text(xmin * 0.9, median_intensity * 0.9,
                f"median rank: {median_pathway_rank} ({median_pathway_rank / len(rank_data) * 100 :.1f}%) "
                f"with intensity: {median_intensity:.2E}",  # TODO case for log and non log
                verticalalignment="top", horizontalalignment="left")

    if "Log_2" not in intensity_label:
        ax.set_yscale("log")
    ax.set_xlabel("Protein rank")
    ax.set_ylabel(f"{full_name} mean")

    fig.legend(bbox_to_anchor=(1.02, 0.5), loc="center left")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("{sample1}_vs_{sample2}")
def save_experiment_comparison_results(
        protein_intensities_sample1: pd.Series, protein_intensities_sample2: pd.Series,
        exclusive_sample1: pd.Series, exclusive_sample2: pd.Series, sample1: str, sample2: str,
        intensity_label: str = "Intensity", show_suptitle: bool = False,
        plot: Optional[Tuple[plt.Figure, plt.Axes]] = None,  **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    # calculate r
    try:
        r = stats.pearsonr(protein_intensities_sample1, protein_intensities_sample2)
    except ValueError:
        warnings.warn(f"Could not calculate pearson r for {sample1} vs {sample2}")
        r = (np.nan,)

    if plot is not None:
        fig, ax = plot
    else:
        plt.close("all")
        fig, ax = plt.subplots(1, 1, figsize=(7, 7))

    exp = r"$r^{2}$"
    ax.scatter(protein_intensities_sample1, protein_intensities_sample2, s=8, alpha=0.6, marker=".",
               label=f"{sample1} vs {sample2}, {exp}: {r[0] ** 2:.4f}")
    ax.scatter(exclusive_sample1, [np.min(protein_intensities_sample2) * 0.95] * exclusive_sample1.shape[0],
               s=8, alpha=0.6, marker=".", label=f"exclusive for {sample1}")
    ax.scatter([np.min(protein_intensities_sample1) * 0.95] * exclusive_sample2.shape[0], exclusive_sample2,
               s=8, alpha=0.6, marker=".", label=f"exclusive for {sample2}")

    ax.set_xlabel(sample1)
    ax.set_ylabel(sample2)
    fig.legend(frameon=False)
    if "Log_2" not in intensity_label:
        ax.set_xscale("log")
        ax.set_yscale("log")

    xmin, xmax = ax.get_xbound()
    ymin, ymax = ax.get_ybound()
    ax.set_xlim(min(xmin, ymin), max(xmax, ymax))
    ax.set_ylim(min(xmin, ymin), max(xmax, ymax))

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("go_analysis")
def save_go_analysis_results(
        heights, test_results, go_analysis_gene_names, intensity_label="Intensity", **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    # TODO also create table
    # TODO move labels to bars, remove legend
    plt.close("all")
    fig, ax = plt.subplots(1, 1, figsize=(7, int(len(heights) * len(go_analysis_gene_names) / 3)))

    bar_width = 0.25
    for i, experiment in enumerate(heights):
        y_pos = np.array([(i + (len(heights) + 1) * x) * bar_width for x in range(len(go_analysis_gene_names))])
        ax.barh(y_pos, heights[experiment], height=bar_width, edgecolor='white', label=experiment)
        for x, y, text in zip(heights[experiment], y_pos, test_results[experiment]):
            if text > 0.05:
                continue
            text = f"{text:.4f}" if text > 0.0005 else "< 0.0005"
            ax.text(x, y, f" p: {text}", verticalalignment="center", fontsize="8")

    ax.set_ylabel('compartiment')
    ax.set_xlabel('number of proteins')
    # add yticks on the middle of the group bars
    ax.set_yticks([(len(heights) / 2 + (len(heights) + 1) * x) * bar_width
                   for x in range(len(go_analysis_gene_names))])
    # replace the y ticks with the compartiment names
    ax.set_yticklabels([x for x in go_analysis_gene_names])
    plt.legend()

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("pathway_timecourse_{pathway}")
def save_pathway_timecourse_results():
    """plt.close("all")
    n_rows, n_cols = get_number_rows_cols_for_fig(found_proteins)
    fig, axarr = plt.subplots(n_rows, n_cols, figsize=(n_cols * int(max_time / 5), 4 * n_rows))
    if show_suptitle:
        fig.suptitle(pathway)

    protein_minimum = self.all_intensities_dict[df_to_use].max().max()
    protein_maximum = self.all_intensities_dict[df_to_use].min().min()
    for protein, (pos, ax) in zip(found_proteins, np.ndenumerate(axarr)):
        ax.set_title(protein)
        ax.set_xlabel(f"Age [weeks]")
        ax.set_ylabel(intensity_label)
        for idx, experiment in enumerate(level_keys):
            protein_intensities = self.all_tree_dict[df_to_use][experiment].aggregate(None, index=protein)
            mask = protein_intensities > 0
            protein_minimum = min(protein_minimum, protein_intensities[mask].min())
            protein_maximum = max(protein_maximum, protein_intensities[mask].max())
            ax.scatter([x_values[experiment]] * sum(mask), protein_intensities[mask],
                       label=f"{groups[experiment]}", color=group_colors[groups[experiment]])
    # adjust labels based on overall min and max of the pathway
    for protein, (pos, ax) in zip(found_proteins, np.ndenumerate(axarr)):
        ax.set_ylim(bottom=protein_minimum * 0.99, top=protein_maximum * 1.01)
        ax.set_xlim(left=0, right=max_time + 1)
    handles, labels = next(np.ndenumerate(axarr))[1].get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(1.04, 0.5), loc="center left")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    return fig, axarr"""


@save_plot("venn_bar_{ex}")
@save_venn_to_txt({"named_sets": "set_bar"})
def save_bar_venn(
        named_sets: Dict[str, set], ex: str, show_suptitle: bool = True, **kwargs
) -> Optional[Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]]:
    plt.close("all")
    if len(named_sets) > 6:
        warnings.warn(f"Skipping bar-venn for {ex} because it has more than 6 experiments")
        return

    # create a mapping from name to a y coordinate
    y_mappings = {name: i for i, name in enumerate(named_sets)}
    # get all the heights and other info required for the plot
    heights = []
    x = []
    ys = []
    for i, (intersected, unioned, result) in enumerate(venn_names(named_sets)):
        heights.append(len(result))
        x.append(i)
        ys.append([y_mappings[x] for x in intersected])

    # initial figure setup
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(1 * len(heights), 7))
    if show_suptitle:
        fig.suptitle(ex, fontsize=20)
    # create the bar plot
    ax1.bar(x, heights, color="skyblue")
    # add text to the bar plot
    for x_level, height in zip(x, heights):
        ax1.text(x_level, max(heights) / 2, height, verticalalignment='center', horizontalalignment='center')
    ax1.set_ylabel("Number of proteins")

    # create the line plots
    for x_level, y in zip(x, ys):
        # we just want to draw a straight line every time so we repeat x as often as needed
        ax2.plot([x_level] * len(y), y, linestyle="-", color="gray", marker=".")
    # replace the yticks with the names of the samples
    ax2.set_yticks([i for i in range(len(y_mappings))])
    ax2.set_yticklabels(y_mappings)
    ax2.set_ylabel("Sample name")
    ax2.set_xlabel("Number of comparison")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, (ax1, ax2)


@save_plot("venn_replicate_{ex}")
@save_venn_to_txt({"named_sets": "set"})
def save_venn(
        named_sets: Dict[str, set], ex: str, show_suptitle: bool = True,
        title_font_size=20, set_label_font_size=16, subset_label_font_size=14, **kwargs
) -> Optional[Tuple[plt.Figure, plt.Axes]]:
    plt.close("all")
    fig, ax = plt.subplots(1, 1, figsize=(14, 7))
    if show_suptitle:
        plt.title(ex, fontsize=title_font_size)

    # create venn diagram based on size of set
    sets = named_sets.values()
    set_names = named_sets.keys()
    if len(sets) < 2:
        warnings.warn(f"Could not create venn diagram for {ex} because it has less than 2 replicates")
        return
    elif len(sets) == 2:
        venn = venn2(subsets=sets, set_labels=set_names, ax=ax)
    elif len(sets) == 3:
        venn = venn3(subsets=sets, set_labels=set_names, ax=ax)
    else:
        warnings.warn(f"Could not create venn diagram for {ex}"
                      f" because it has more than 3 replicates ({len(sets)})")
        return

    # if a figure was created, do some further configuration
    for text in venn.set_labels:
        try:
            text.set_fontsize(set_label_font_size)
        except AttributeError:
            pass
    handles = []
    labels = []
    for text, patch in zip(venn.subset_labels, venn.patches):
        try:
            handles.append(patch)
            labels.append(text.get_text())
            text.set_fontsize(subset_label_font_size)
        except AttributeError:
            pass
    plt.legend(handles, labels, bbox_to_anchor=(1.02, 0.5), loc="center left")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax
